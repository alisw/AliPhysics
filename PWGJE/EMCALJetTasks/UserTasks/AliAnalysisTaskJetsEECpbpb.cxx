//Task for ENC in pbpb collisions
//Code inherited and modified from AliAnalysisTaskLundPlane.cxx authored by Leticia Cunqueiro and Laura Havener
//Authors: Ananya Rai, Laura Havener
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetContainer.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliMCParticleContainer.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliTLorentzVector.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisDataSlot.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TMath.h>
#include <TGrid.h>
#include "AliAODEvent.h"
#include "AliEmcalContainerUtils.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliVParticle.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskRho.h"
#include "AliVHeader.h"
#include "AliAODVertex.h"
#include "AliFJWrapper.h"
#include "AliAODPid.h"
#include "AliFJWrapper.h"
#include "AliAnalysisTaskJetsEECpbpb.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetsEECpbpb)

//________________________________________________________________________
AliAnalysisTaskJetsEECpbpb::AliAnalysisTaskJetsEECpbpb(): AliAnalysisTaskEmcalJet("AliAnalysisTaskJetsEECpbpb", kTRUE),
fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(5), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.9), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fMatchR(0.3), fjetMinPtSub(20), fjetMaxPtSub(200), fjetMinArea(0), fEmbTuthJetPtMin(0), fStoreTrig(kFALSE), fConeR(0.4), fpTcorr(0), fpaircut(0), 
fpairfastsim(0), fUnfolding(0),fMatchJetTrack(1),fMissJetTrack(1),fFakeJetTrack(1), fMaxPtTrack(100), fJetPtMin(0),
fDoEmbedding(kFALSE), fIsEmbeddedEvent(kFALSE),fDoPerpCone(kFALSE),fDoRandCone(kFALSE), fCout(kTRUE), fisGoodIncEvent(0),
fDoPartLevelMatching(0),fDoDetLevelMatching(0),fTruthMinLabel(0),fTruthMaxLabel(0),fSaveMCInformation(0),fJetMatchingSharedPtFraction(0),fjetRhoVal(0),fjetRecoRhoVal(0),
fjetEmbRhoVal(0),fjetGenRhoVal(0),fhasAcceptedEMCjet(0),fJetContainer(0),
fbgJetContainer(0),fRecoJetContainer(0),fEmbJetContainer(0),fGenJetContainer(0),fGeneratorLevelName(), fDetectorLevelName(), fGeneratorLevel(0), fDetectorLevel(0), 
fJet_truCont(0), fAOD(0),jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), E3C_hist(0), E3C_pt_hist(0),  EEC_det_pt_hist_3d(0), 
EEC_tru_pt_hist_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), 
N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), 
E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), pt_tru(0), pt_tru_match(0), pt_det(0), 
pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), 
pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0), qpt_tru(0), qpt_det(0), track_pt_tru(0), 
track_pt_det(0), track_pt_matched(0), R_match_eec(0), wt_match_eec(0), R_match_e3c(0), 
wt_match_e3c(0), qpt_tru1(0), qpt_tru2(0),pt_tru1(0), pt_tru2(0), eec_Mm(0), eec_mm(0), 
e3c_MMm(0), e3c_Mmm(0), e3c_mmm(0), eec_Mf(0), eec_ff(0), e3c_MMf(0), e3c_Mff(0), e3c_fff(0), 
eec_matched_det(0), eec_matched_tru(0), e3c_matched_det(0), e3c_matched_tru(0), wt_res_eec(0), 
wt_res_e3c(0), R_res_eec(0), R_res_e3c(0), wtnojet_match_eec(0), wtnojet_match_e3c(0), 
wtnojet_res_eec(0), wtnojet_res_e3c(0), R_match_eec_tru(0), wt_match_eec_tru(0), R_match_e3c_tru(0), 
wt_match_e3c_tru(0), wt_res_eec_tru(0), wt_res_e3c_tru(0), R_res_eec_tru(0), R_res_e3c_tru(0), 
wtnojet_match_eec_tru(0), wtnojet_match_e3c_tru(0), wtnojet_res_eec_tru(0), wtnojet_res_e3c_tru(0), 
constituentId(0),R_res_eec_tru_debug(0),track_pt_res_debug(0),track_eta_debug(0), track_rap_debug(0), 
track_phi_debug(0), track_R_debug(0), track_R_debug_rap(0), track_eta_res_debug(0), track_rap_res_debug(0), 
track_phi_res_debug(0), track_R_res_debug(0), track_R_rap_res_debug(0), R_match_eec_tru_rap(0), 
track_pt_response_debug(0), track_pt_wt_res_debug(0), track_pt_wt_response_debug(0), jetpt_res_w_R(0), jetpt_res_w_wt(0),
NumJetEvent(0), h_JMB_dat(0), h_MB1_dat(0), h_MB1MB2_dat(0), h3Jet_deltaR_JMB_dat(0), h3Jet_deltaR_MB1_dat(0),h3Jet_deltaR_MB1MB2_dat(0), 
h_MJ(0), h_MJ0(0), h_MJ1(0), h_MJ2(0), h_MJ_m(0), h_MJ0_m(0), h_MJ1_m(0), h_MJ2_m(0), h_MJ_tru_m(0), h_MJ0_tru_m(0),
h_MJ1_tru_m(0), h_MJ2_tru_m(0), h_MJ_um(0), h_MJ0_um(0), h_MJ1_um(0), h_MJ2_um(0), h_MJ_c(0), h_MJ0_c(0), h_MJ1_c(0),
h_MJ2_c(0), h_MJ_c_m(0), h_MJ0_c_m(0), h_MJ1_c_m(0), h_MJ2_c_m(0), h_MJ_tru_c_m(0), h_MJ0_tru_c_m(0), h_MJ1_tru_c_m(0),
h_MJ2_tru_c_m(0), h_MJ_c_um(0), h_MJ0_c_um(0), h_MJ1_c_um(0), h_MJ2_c_um(0), h_MB1(0), h_MB1_m(0), h_MB1_tru_m(0),
h_MB1_um(0), h_MB1_c(0), h_MB1_c_m(0), h_MB1_tru_c_m(0), h_MB1_c_um(0), h_JMB(0), h_JMB_m(0),
h_JMB_tru_m(0), h_JMB_um(0), h_JMB_c(0), h_JMB_c_m(0), h_JMB_tru_c_m(0), h_JMB_c_um(0), 
h_SMB(0), h_SMB_m(0), h_SMB_tru_m(0), h_SMB_um(0), h_SMB_c(0), h_SMB_c_m(0), h_SMB_tru_c_m(0),
h_SMB_c_um(0), h_BMB(0), h_BMB_m(0), h_BMB_tru_m(0), h_BMB_um(0), h_BMB_c(0), h_BMB_c_m(0),
h_BMB_tru_c_m(0), h_BMB_c_um(0), h_MB1MB2(0), h_MB1MB2_m(0), h_MB1MB2_tru_m(0), h_MB1MB2_um(0),
h_MB1MB2_c(0), h_MB1MB2_c_m(0), h_MB1MB2_tru_c_m(0), h_MB1MB2_c_um(0), h_MJ_e3c(0), h_MJ0_e3c(0),
h_MJ1_e3c(0), h_MJ2_e3c(0), h_MJ3_e3c(0), h_MJ_e3c_tru(0), h_MJ0_e3c_tru(0), h_MJ1_e3c_tru(0),
h_MJ2_e3c_tru(0), h_MJ3_e3c_tru(0), h_MB1MB1MB1(0), h_MB1MB1MB1_dat(0), h_MB1MB1MB1_tru(0), h_JJMB(0),
h_JJMB_dat(0), h_JJMB_tru(0), h_JMBMB(0), h_JMBMB_dat(0), h_JMBMB_tru(0), h_MB1MB1MB2(0), h_MB1MB1MB2_dat(0),
h_MB1MB1MB2_tru(0), h_MB1MB2MB2(0), h_MB1MB2MB2_dat(0), h_MB1MB2MB2_tru(0), h_JMB1MB2(0), h_JMB1MB2_dat(0),
h_JMB1MB2_tru(0), h_MB1MB2MB3(0), h_MB1MB2MB3_dat(0), h_MB1MB2MB3_tru(0), h_BMBMB(0), h_SMBMB(0), h_BBMB(0),
h_BMB1MB2(0), h_SMB1MB2(0), h_SBMB(0), h_SSMB(0), h_BBMB_tru(0), h_SBMB_tru(0), h_SSMB_tru(0), h_BMBMB_tru(0),
h_SMBMB_tru(0), h_BMB1MB2_tru(0), h_SMB1MB2_tru(0), h3Jet_deltaR_MJ(0), h3Jet_deltaR_MJ0(0), h3Jet_deltaR_MJ1(0),
h3Jet_deltaR_MJ2(0), h3Jet_deltaR_MJ_m(0), h3Jet_deltaR_MJ0_m(0), h3Jet_deltaR_MJ1_m(0), h3Jet_deltaR_MJ2_m(0),
h3Jet_deltaR_MJ_tru_m(0), h3Jet_deltaR_MJ0_tru_m(0), h3Jet_deltaR_MJ1_tru_m(0), h3Jet_deltaR_MJ2_tru_m(0),
h3Jet_deltaR_MJ_minpt(0), h3Jet_deltaR_MJ0_minpt(0), h3Jet_deltaR_MJ1_minpt(0), h3Jet_deltaR_MJ2_minpt(0),
h3Jet_deltaR_MJ_tru_minpt(0), h3Jet_deltaR_MJ0_tru_minpt(0), h3Jet_deltaR_MJ1_tru_minpt(0), h3Jet_deltaR_MJ2_tru_minpt(0),
h3Jet_deltaR_MJ_um(0), h3Jet_deltaR_MJ0_um(0), h3Jet_deltaR_MJ1_um(0), h3Jet_deltaR_MJ2_um(0), h3Jet_deltaR_MJ_c(0),
h3Jet_deltaR_MJ0_c(0), h3Jet_deltaR_MJ1_c(0), h3Jet_deltaR_MJ2_c(0), h3Jet_deltaR_MJ_c_m(0), h3Jet_deltaR_MJ0_c_m(0),
h3Jet_deltaR_MJ1_c_m(0), h3Jet_deltaR_MJ2_c_m(0), h3Jet_deltaR_MJ_tru_c_m(0), h3Jet_deltaR_MJ0_tru_c_m(0), h3Jet_deltaR_MJ1_tru_c_m(0),
h3Jet_deltaR_MJ2_tru_c_m(0), h3Jet_deltaR_MB1(0), h3Jet_deltaR_MB1_m(0), h3Jet_deltaR_MB1_tru_m(0), h3Jet_deltaR_MB1_minpt(0),
h3Jet_deltaR_MB1_tru_minpt(0),h3Jet_deltaR_MB1_um(0), h3Jet_deltaR_MB1_c(0), h3Jet_deltaR_MB1_c_m(0), h3Jet_deltaR_MB1_tru_c_m(0),
h3Jet_deltaR_MB1_c_um(0), h3Jet_deltaR_JMB(0), h3Jet_deltaR_JMB_m(0), h3Jet_deltaR_JMB_tru_m(0), h3Jet_deltaR_JMB_c(0), h3Jet_deltaR_JMB_c_m(0),
h3Jet_deltaR_JMB_c_um(0), h3Jet_deltaR_SMB(0), h3Jet_deltaR_SMB_m(0), h3Jet_deltaR_SMB_tru_m(0), h3Jet_deltaR_SMB_c(0), h3Jet_deltaR_SMB_c_m(0),
h3Jet_deltaR_SMB_c_um(0), OptUn_eec(0), OptUn_e3c(0), hJet_num_um(0), hJet_num_m(0), h_jetpt_m(0), h_jetpt_um(0), h_jetpt_m_tru(0),
h3_MB1MB1MB1_dat(0),h3_JJMB_dat(0),h3_JMB1MB2_dat(0),h3_JMBMB_dat(0),h3_MB1MB1MB2_dat(0),h3_MB1MB2MB2_dat(0),h3_MB1MB2MB3_dat(0),
hJet_deltaR_MJ_e3c(0),hJet_deltaR_MJ0_e3c(0),hJet_deltaR_MJ1_e3c(0),hJet_deltaR_MJ2_e3c(0),hJet_deltaR_MJ3_e3c(0),
hJet_deltaR_MJ_e3c_m(0),hJet_deltaR_MJ0_e3c_m(0),hJet_deltaR_MJ1_e3c_m(0),hJet_deltaR_MJ2_e3c_m(0),hJet_deltaR_MJ3_e3c_m(0),h_MB1MB1MB1_m(0),
h_JJMB_m(0),h_JMBMB_m(0),h_MB1MB1MB2_m(0),h_MB1MB2MB2_m(0),h_JMB1MB2_m(0),h_MB1MB2MB3_m(0),h_BMBMB_m(0),h_SMBMB_m(0),h_BMB1MB2_m(0),
h_SMB1MB2_m(0),h_BBMB_m(0),h_SBMB_m(0),h_SSMB_m(0),
hJet_deltaR_MJ_e3c_tru_m(0),hJet_deltaR_MJ0_e3c_tru_m(0),hJet_deltaR_MJ1_e3c_tru_m(0),hJet_deltaR_MJ2_e3c_tru_m(0),hJet_deltaR_MJ3_e3c_tru_m(0),
h_MB1MB1MB1_tru_m(0),h_JJMB_tru_m(0),h_JMBMB_tru_m(0),h_MB1MB1MB2_tru_m(0),
h_MB1MB2MB2_tru_m(0),h_JMB1MB2_tru_m(0),h_MB1MB2MB3_tru_m(0),h_BMBMB_tru_m(0),h_SMBMB_tru_m(0),h_BMB1MB2_tru_m(0),
h_BBMB_tru_m(0),h_SBMB_tru_m(0),h_SSMB_tru_m(0),h_SMB1MB2_tru_m(0),
hJet_deltaR_MJ_e3c_um(0),hJet_deltaR_MJ0_e3c_um(0),hJet_deltaR_MJ1_e3c_um(0),
hJet_deltaR_MJ2_e3c_um(0),hJet_deltaR_MJ3_e3c_um(0),h_MB1MB1MB1_um(0),h_JJMB_um(0),h_JMBMB_um(0),h_MB1MB1MB2_um(0),h_MB1MB2MB2_um(0),
h_JMB1MB2_um(0),h_MB1MB2MB3_um(0),h_BMBMB_um(0),h_SMBMB_um(0),h_BMB1MB2_um(0),h_SMB1MB2_um(0),h_BBMB_um(0),h_SBMB_um(0),h_SSMB_um(0),
h3Jet_deltaR_MJ_e3c(0),h3Jet_deltaR_MJ0_e3c(0),h3Jet_deltaR_MJ1_e3c(0),h3Jet_deltaR_MJ2_e3c(0),h3Jet_deltaR_MJ3_e3c(0),h3_MB1MB1MB1(0),
h3_JJMB(0),h3_JMBMB(0),h3_MB1MB1MB2(0),h3_MB1MB2MB2(0),h3_JMB1MB2(0),h3_MB1MB2MB3(0),h3_BMBMB(0),h3_SMBMB(0),h3_BMB1MB2(0),
h3_SMB1MB2(0),h3_BBMB(0),h3_SBMB(0),h3_SSMB(0),h3Jet_deltaR_MJ_e3c_m(0),h3Jet_deltaR_MJ1_e3c_m(0),h3Jet_deltaR_MJ2_e3c_m(0),
h3Jet_deltaR_MJ3_e3c_m(0),h3_MB1MB1MB1_m(0),h3_JJMB_m(0),h3_JMBMB_m(0),h3_MB1MB1MB2_m(0),h3_MB1MB2MB2_m(0),h3_JMB1MB2_m(0),
h3_MB1MB2MB3_m(0),h3_BMBMB_m(0),h3_SMBMB_m(0),h3_BMB1MB2_m(0),h3_SMB1MB2_m(0),h3_BBMB_m(0),h3_SBMB_m(0),h3_SSMB_m(0),
h3Jet_deltaR_MJ_e3c_tru_m(0),h3Jet_deltaR_MJ1_e3c_tru_m(0),h3Jet_deltaR_MJ2_e3c_tru_m(0),h3Jet_deltaR_MJ3_e3c_tru_m(0),
h3_MB1MB1MB1_tru_m(0),h3_JJMB_tru_m(0),h3_JMBMB_tru_m(0),h3_MB1MB1MB2_tru_m(0),h3_MB1MB2MB2_tru_m(0),h3_JMB1MB2_tru_m(0),
h3_MB1MB2MB3_tru_m(0),h3_BMBMB_tru_m(0),h3_SMBMB_tru_m(0),h3_BMB1MB2_tru_m(0),h3_SMB1MB2_tru_m(0),h3_BBMB_tru_m(0),h3_SBMB_tru_m(0),
h3_SSMB_tru_m(0),h3Jet_deltaR_MJ3_e3c_um(0),h3_MB1MB1MB1_um(0),h3_JJMB_um(0),h3_JMBMB_um(0),h3_MB1MB1MB2_um(0),h3_MB1MB2MB2_um(0),
h3_JMB1MB2_um(0),h3_MB1MB2MB3_um(0),h3_BMBMB_um(0),h3_SMBMB_um(0),h3_BMB1MB2_um(0),h3_SMB1MB2_um(0),h3_BBMB_um(0),h3_SBMB_um(0),
h3_SSMB_um(0),fTreeMatchTracks(0), fTreeData(0), fMCParticleArrayName("mcparticles"), fMCParticleArray(0),
fRandom(0x0), ifeec(1), ife3c(0), ifMinPtHist(0), ifcFactorHist(0), fAddEventCuts(0), fHighPtTrackCutEvent(0),fDeltaAxisShift(0.6),
h_dpt(0),delta_pt_cone(0),delta_pt_coneBias1(0),delta_pt_coneBias2(0),delta_pt_ENCcone(0),delta_pt_coneBias3(0),delta_pt_coneBias4(0)
{
  if(fCout){
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  std::cout << " Info::anrai:***************** CONSTRUCTOR CALLED: AliAnalysisTaskJetsEECpbpb *******************"<< std::endl;
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEECpbpb::AliAnalysisTaskJetsEECpbpb(const char *name): AliAnalysisTaskEmcalJet(name, kTRUE),
fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(5), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.9), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fMatchR(0.3), fjetMinPtSub(20), fjetMaxPtSub(200), fjetMinArea(0), fEmbTuthJetPtMin(0), fStoreTrig(kFALSE), fConeR(0.4), fpTcorr(0), fpaircut(0), 
fpairfastsim(0), fUnfolding(0),fMatchJetTrack(1),fMissJetTrack(1),fFakeJetTrack(1), fMaxPtTrack(100), fJetPtMin(0),
fDoEmbedding(kFALSE), fIsEmbeddedEvent(kFALSE),fDoPerpCone(kFALSE),fDoRandCone(kFALSE), fCout(kTRUE), fisGoodIncEvent(0),
fDoPartLevelMatching(0),fDoDetLevelMatching(0),fTruthMinLabel(0),fTruthMaxLabel(0),fSaveMCInformation(0),fJetMatchingSharedPtFraction(0),fjetRhoVal(0),fjetRecoRhoVal(0),
fjetEmbRhoVal(0),fjetGenRhoVal(0), fhasAcceptedEMCjet(0),fJetContainer(0),
fbgJetContainer(0),fRecoJetContainer(0),fEmbJetContainer(0),fGenJetContainer(0),fGeneratorLevelName(), fDetectorLevelName(), fGeneratorLevel(0), fDetectorLevel(0), 
fJet_truCont(0), fAOD(0),jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), E3C_hist(0), E3C_pt_hist(0),  EEC_det_pt_hist_3d(0), 
EEC_tru_pt_hist_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), 
N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), 
E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), pt_tru(0), pt_tru_match(0), pt_det(0), 
pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), 
pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0), qpt_tru(0), qpt_det(0), track_pt_tru(0), 
track_pt_det(0), track_pt_matched(0), R_match_eec(0), wt_match_eec(0), R_match_e3c(0), 
wt_match_e3c(0), qpt_tru1(0), qpt_tru2(0),pt_tru1(0), pt_tru2(0), eec_Mm(0), eec_mm(0), 
e3c_MMm(0), e3c_Mmm(0), e3c_mmm(0), eec_Mf(0), eec_ff(0), e3c_MMf(0), e3c_Mff(0), e3c_fff(0), 
eec_matched_det(0), eec_matched_tru(0), e3c_matched_det(0), e3c_matched_tru(0), wt_res_eec(0), 
wt_res_e3c(0), R_res_eec(0), R_res_e3c(0), wtnojet_match_eec(0), wtnojet_match_e3c(0), 
wtnojet_res_eec(0), wtnojet_res_e3c(0), R_match_eec_tru(0), wt_match_eec_tru(0), R_match_e3c_tru(0), 
wt_match_e3c_tru(0), wt_res_eec_tru(0), wt_res_e3c_tru(0), R_res_eec_tru(0), R_res_e3c_tru(0), 
wtnojet_match_eec_tru(0), wtnojet_match_e3c_tru(0), wtnojet_res_eec_tru(0), wtnojet_res_e3c_tru(0), 
constituentId(0),R_res_eec_tru_debug(0),track_pt_res_debug(0),track_eta_debug(0), track_rap_debug(0), 
track_phi_debug(0), track_R_debug(0), track_R_debug_rap(0), track_eta_res_debug(0), track_rap_res_debug(0), 
track_phi_res_debug(0), track_R_res_debug(0), track_R_rap_res_debug(0), R_match_eec_tru_rap(0), 
track_pt_response_debug(0), track_pt_wt_res_debug(0), track_pt_wt_response_debug(0), jetpt_res_w_R(0), jetpt_res_w_wt(0),
NumJetEvent(0), h_JMB_dat(0), h_MB1_dat(0), h_MB1MB2_dat(0), h3Jet_deltaR_JMB_dat(0), h3Jet_deltaR_MB1_dat(0),h3Jet_deltaR_MB1MB2_dat(0), 
h_MJ(0), h_MJ0(0), h_MJ1(0), h_MJ2(0), h_MJ_m(0), h_MJ0_m(0), h_MJ1_m(0), h_MJ2_m(0), h_MJ_tru_m(0), h_MJ0_tru_m(0),
h_MJ1_tru_m(0), h_MJ2_tru_m(0), h_MJ_um(0), h_MJ0_um(0), h_MJ1_um(0), h_MJ2_um(0), h_MJ_c(0), h_MJ0_c(0), h_MJ1_c(0),
h_MJ2_c(0), h_MJ_c_m(0), h_MJ0_c_m(0), h_MJ1_c_m(0), h_MJ2_c_m(0), h_MJ_tru_c_m(0), h_MJ0_tru_c_m(0), h_MJ1_tru_c_m(0),
h_MJ2_tru_c_m(0), h_MJ_c_um(0), h_MJ0_c_um(0), h_MJ1_c_um(0), h_MJ2_c_um(0), h_MB1(0), h_MB1_m(0), h_MB1_tru_m(0),
h_MB1_um(0), h_MB1_c(0), h_MB1_c_m(0), h_MB1_tru_c_m(0), h_MB1_c_um(0), h_JMB(0), h_JMB_m(0),
h_JMB_tru_m(0), h_JMB_um(0), h_JMB_c(0), h_JMB_c_m(0), h_JMB_tru_c_m(0), h_JMB_c_um(0), 
h_SMB(0), h_SMB_m(0), h_SMB_tru_m(0), h_SMB_um(0), h_SMB_c(0), h_SMB_c_m(0), h_SMB_tru_c_m(0),
h_SMB_c_um(0), h_BMB(0), h_BMB_m(0), h_BMB_tru_m(0), h_BMB_um(0), h_BMB_c(0), h_BMB_c_m(0),
h_BMB_tru_c_m(0), h_BMB_c_um(0), h_MB1MB2(0), h_MB1MB2_m(0), h_MB1MB2_tru_m(0), h_MB1MB2_um(0),
h_MB1MB2_c(0), h_MB1MB2_c_m(0), h_MB1MB2_tru_c_m(0), h_MB1MB2_c_um(0), h_MJ_e3c(0), h_MJ0_e3c(0),
h_MJ1_e3c(0), h_MJ2_e3c(0), h_MJ3_e3c(0), h_MJ_e3c_tru(0), h_MJ0_e3c_tru(0), h_MJ1_e3c_tru(0),
h_MJ2_e3c_tru(0), h_MJ3_e3c_tru(0), h_MB1MB1MB1(0), h_MB1MB1MB1_dat(0), h_MB1MB1MB1_tru(0), h_JJMB(0),
h_JJMB_dat(0), h_JJMB_tru(0), h_JMBMB(0), h_JMBMB_dat(0), h_JMBMB_tru(0), h_MB1MB1MB2(0), h_MB1MB1MB2_dat(0),
h_MB1MB1MB2_tru(0), h_MB1MB2MB2(0), h_MB1MB2MB2_dat(0), h_MB1MB2MB2_tru(0), h_JMB1MB2(0), h_JMB1MB2_dat(0),
h_JMB1MB2_tru(0), h_MB1MB2MB3(0), h_MB1MB2MB3_dat(0), h_MB1MB2MB3_tru(0), h_BMBMB(0), h_SMBMB(0), h_BBMB(0),
h_BMB1MB2(0), h_SMB1MB2(0), h_SBMB(0), h_SSMB(0), h_BBMB_tru(0), h_SBMB_tru(0), h_SSMB_tru(0), h_BMBMB_tru(0),
h_SMBMB_tru(0), h_BMB1MB2_tru(0), h_SMB1MB2_tru(0), h3Jet_deltaR_MJ(0), h3Jet_deltaR_MJ0(0), h3Jet_deltaR_MJ1(0),
h3Jet_deltaR_MJ2(0), h3Jet_deltaR_MJ_m(0), h3Jet_deltaR_MJ0_m(0), h3Jet_deltaR_MJ1_m(0), h3Jet_deltaR_MJ2_m(0),
h3Jet_deltaR_MJ_tru_m(0), h3Jet_deltaR_MJ0_tru_m(0), h3Jet_deltaR_MJ1_tru_m(0), h3Jet_deltaR_MJ2_tru_m(0),
h3Jet_deltaR_MJ_minpt(0), h3Jet_deltaR_MJ0_minpt(0), h3Jet_deltaR_MJ1_minpt(0), h3Jet_deltaR_MJ2_minpt(0),
h3Jet_deltaR_MJ_tru_minpt(0), h3Jet_deltaR_MJ0_tru_minpt(0), h3Jet_deltaR_MJ1_tru_minpt(0), h3Jet_deltaR_MJ2_tru_minpt(0),
h3Jet_deltaR_MJ_um(0), h3Jet_deltaR_MJ0_um(0), h3Jet_deltaR_MJ1_um(0), h3Jet_deltaR_MJ2_um(0), h3Jet_deltaR_MJ_c(0),
h3Jet_deltaR_MJ0_c(0), h3Jet_deltaR_MJ1_c(0), h3Jet_deltaR_MJ2_c(0), h3Jet_deltaR_MJ_c_m(0), h3Jet_deltaR_MJ0_c_m(0),
h3Jet_deltaR_MJ1_c_m(0), h3Jet_deltaR_MJ2_c_m(0), h3Jet_deltaR_MJ_tru_c_m(0), h3Jet_deltaR_MJ0_tru_c_m(0), h3Jet_deltaR_MJ1_tru_c_m(0),
h3Jet_deltaR_MJ2_tru_c_m(0), h3Jet_deltaR_MB1(0), h3Jet_deltaR_MB1_m(0), h3Jet_deltaR_MB1_tru_m(0), h3Jet_deltaR_MB1_minpt(0),
h3Jet_deltaR_MB1_tru_minpt(0),h3Jet_deltaR_MB1_um(0), h3Jet_deltaR_MB1_c(0), h3Jet_deltaR_MB1_c_m(0), h3Jet_deltaR_MB1_tru_c_m(0),
h3Jet_deltaR_MB1_c_um(0), h3Jet_deltaR_JMB(0), h3Jet_deltaR_JMB_m(0), h3Jet_deltaR_JMB_tru_m(0), h3Jet_deltaR_JMB_c(0), h3Jet_deltaR_JMB_c_m(0),
h3Jet_deltaR_JMB_c_um(0), h3Jet_deltaR_SMB(0), h3Jet_deltaR_SMB_m(0), h3Jet_deltaR_SMB_tru_m(0), h3Jet_deltaR_SMB_c(0), h3Jet_deltaR_SMB_c_m(0),
h3Jet_deltaR_SMB_c_um(0), OptUn_eec(0), OptUn_e3c(0), hJet_num_um(0), hJet_num_m(0), h_jetpt_m(0), h_jetpt_um(0), h_jetpt_m_tru(0),
h3_MB1MB1MB1_dat(0),h3_JJMB_dat(0),h3_JMB1MB2_dat(0),h3_JMBMB_dat(0),h3_MB1MB1MB2_dat(0),h3_MB1MB2MB2_dat(0),h3_MB1MB2MB3_dat(0),
hJet_deltaR_MJ_e3c(0),hJet_deltaR_MJ0_e3c(0),hJet_deltaR_MJ1_e3c(0),hJet_deltaR_MJ2_e3c(0),hJet_deltaR_MJ3_e3c(0),
hJet_deltaR_MJ_e3c_m(0),hJet_deltaR_MJ0_e3c_m(0),hJet_deltaR_MJ1_e3c_m(0),hJet_deltaR_MJ2_e3c_m(0),hJet_deltaR_MJ3_e3c_m(0),h_MB1MB1MB1_m(0),
h_JJMB_m(0),h_JMBMB_m(0),h_MB1MB1MB2_m(0),h_MB1MB2MB2_m(0),h_JMB1MB2_m(0),h_MB1MB2MB3_m(0),h_BMBMB_m(0),h_SMBMB_m(0),h_BMB1MB2_m(0),
h_SMB1MB2_m(0),h_BBMB_m(0),h_SBMB_m(0),h_SSMB_m(0),
hJet_deltaR_MJ_e3c_tru_m(0),hJet_deltaR_MJ0_e3c_tru_m(0),hJet_deltaR_MJ1_e3c_tru_m(0),hJet_deltaR_MJ2_e3c_tru_m(0),hJet_deltaR_MJ3_e3c_tru_m(0),
h_MB1MB1MB1_tru_m(0),h_JJMB_tru_m(0),h_JMBMB_tru_m(0),h_MB1MB1MB2_tru_m(0),
h_MB1MB2MB2_tru_m(0),h_JMB1MB2_tru_m(0),h_MB1MB2MB3_tru_m(0),h_BMBMB_tru_m(0),h_SMBMB_tru_m(0),h_BMB1MB2_tru_m(0),
h_BBMB_tru_m(0),h_SBMB_tru_m(0),h_SSMB_tru_m(0),h_SMB1MB2_tru_m(0),
hJet_deltaR_MJ_e3c_um(0),hJet_deltaR_MJ0_e3c_um(0),hJet_deltaR_MJ1_e3c_um(0),
hJet_deltaR_MJ2_e3c_um(0),hJet_deltaR_MJ3_e3c_um(0),h_MB1MB1MB1_um(0),h_JJMB_um(0),h_JMBMB_um(0),h_MB1MB1MB2_um(0),h_MB1MB2MB2_um(0),
h_JMB1MB2_um(0),h_MB1MB2MB3_um(0),h_BMBMB_um(0),h_SMBMB_um(0),h_BMB1MB2_um(0),h_SMB1MB2_um(0),h_BBMB_um(0),h_SBMB_um(0),h_SSMB_um(0),
h3Jet_deltaR_MJ_e3c(0),h3Jet_deltaR_MJ0_e3c(0),h3Jet_deltaR_MJ1_e3c(0),h3Jet_deltaR_MJ2_e3c(0),h3Jet_deltaR_MJ3_e3c(0),h3_MB1MB1MB1(0),
h3_JJMB(0),h3_JMBMB(0),h3_MB1MB1MB2(0),h3_MB1MB2MB2(0),h3_JMB1MB2(0),h3_MB1MB2MB3(0),h3_BMBMB(0),h3_SMBMB(0),h3_BMB1MB2(0),
h3_SMB1MB2(0),h3_BBMB(0),h3_SBMB(0),h3_SSMB(0),h3Jet_deltaR_MJ_e3c_m(0),h3Jet_deltaR_MJ1_e3c_m(0),h3Jet_deltaR_MJ2_e3c_m(0),
h3Jet_deltaR_MJ3_e3c_m(0),h3_MB1MB1MB1_m(0),h3_JJMB_m(0),h3_JMBMB_m(0),h3_MB1MB1MB2_m(0),h3_MB1MB2MB2_m(0),h3_JMB1MB2_m(0),
h3_MB1MB2MB3_m(0),h3_BMBMB_m(0),h3_SMBMB_m(0),h3_BMB1MB2_m(0),h3_SMB1MB2_m(0),h3_BBMB_m(0),h3_SBMB_m(0),h3_SSMB_m(0),
h3Jet_deltaR_MJ_e3c_tru_m(0),h3Jet_deltaR_MJ1_e3c_tru_m(0),h3Jet_deltaR_MJ2_e3c_tru_m(0),h3Jet_deltaR_MJ3_e3c_tru_m(0),
h3_MB1MB1MB1_tru_m(0),h3_JJMB_tru_m(0),h3_JMBMB_tru_m(0),h3_MB1MB1MB2_tru_m(0),h3_MB1MB2MB2_tru_m(0),h3_JMB1MB2_tru_m(0),
h3_MB1MB2MB3_tru_m(0),h3_BMBMB_tru_m(0),h3_SMBMB_tru_m(0),h3_BMB1MB2_tru_m(0),h3_SMB1MB2_tru_m(0),h3_BBMB_tru_m(0),h3_SBMB_tru_m(0),
h3_SSMB_tru_m(0),h3Jet_deltaR_MJ3_e3c_um(0),h3_MB1MB1MB1_um(0),h3_JJMB_um(0),h3_JMBMB_um(0),h3_MB1MB1MB2_um(0),h3_MB1MB2MB2_um(0),
h3_JMB1MB2_um(0),h3_MB1MB2MB3_um(0),h3_BMBMB_um(0),h3_SMBMB_um(0),h3_BMB1MB2_um(0),h3_SMB1MB2_um(0),h3_BBMB_um(0),h3_SBMB_um(0),
h3_SSMB_um(0),fTreeMatchTracks(0), fTreeData(0), fMCParticleArrayName("mcparticles"), fMCParticleArray(0),
fRandom(0x0), ifeec(1), ife3c(0), ifMinPtHist(0), ifcFactorHist(0), fAddEventCuts(0), fHighPtTrackCutEvent(0),fDeltaAxisShift(0.6),
h_dpt(0),delta_pt_cone(0),delta_pt_coneBias1(0),delta_pt_coneBias2(0),delta_pt_ENCcone(0),delta_pt_coneBias3(0),delta_pt_coneBias4(0)
{
  if(fCout){
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  std::cout << " Info::anrai:***************** CONSTRUCTOR CALLED: AliAnalysisTaskJetsEECpbpb *******************"<< std::endl;
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  std::cout << " Info::anrai:===================================================================================="<< std::endl;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

}


//________________________________________________________________________
AliAnalysisTaskJetsEECpbpb::~AliAnalysisTaskJetsEECpbpb() {
  // Destructor.
 if(fCout){std::cout << " Info::anrai:***************** DESTRUCTOR CALLED: AliAnalysisTaskJetsEECpbpb *******************"<< std::endl;}
 if (fRandom) delete fRandom;
}

//________________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::UserCreateOutputObjects() {

if(fCout)  std::cout << " Info::anrai: ===== In the UserCreateOutputObjects ===== " << std::endl;

// Create user output.
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    TH1::AddDirectory(oldStatus);


  TRandom3* fRandom = new TRandom3(0);

  if(fCout) cout<<"Particle collision array"<<endl;
  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++){
    if(fCout) cout<<"In particle collision array loop"<<endl;
    if (GetParticleContainer(iCont)->GetIsEmbedding()){
      fIsEmbeddedEvent = kTRUE;
    }
   }

  if (fDoEmbedding)
  {
    if(!fIsEmbeddedEvent){
      if(fCout) cout<<"In case of non-embedded array"<<endl;
      fMCParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticleArrayName.Data()));
    }
    else
    {
      if(fCout) cout<<"In case of embedding, the MC particle array needs to be fetched differently, fetching now"<<endl;
      // In case of embedding, the MC particle array needs to be fetched differently
      AliVEvent* event = AliEmcalContainerUtils::GetEvent(InputEvent(), kTRUE);
      fMCParticleArray = dynamic_cast<TClonesArray*>(event->FindListObject(fMCParticleArrayName.Data()));
    }
  }

   if(fCout) cout<<"histos being initialized"<<endl;
    Double_t lower = 0.01;
    Double_t upper = 0.47863;
    Int_t bins = 21;  // 21 bins --> 22 bin edges
    // Compute logarithms of the bounds.
    Double_t from = TMath::Log10(lower);  // from = -2
    Double_t to = TMath::Log10(upper);      // to ≈ -0.320
    // Calculate the width of each bin in log space.
    Double_t width = (to - from) / bins;
    Double_t new_bins[22] = {};
    for (Int_t i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10.0, from + i * width);
    }

    // jet pT bins
    Double_t from_const = 10;
    Double_t to_const = 140;
    Int_t bins_const = 13;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t new_bins_const[14] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
    new_bins_const[i] = (from_const + i * width_const);
    }
    
    //x and Y bins for wt with jet pt included --eec_wt and e3c_wt
    Double_t wt_from = 0;
    Double_t wt_to = 1;
    Int_t wt_bins = 200;
    Double_t wt_width = (wt_to-wt_from)/wt_bins;
    Double_t wt_new_bins[201] = {};
    for (Int_t i = 0; i <= wt_bins; i++)
    {
        wt_new_bins[i] = (wt_from + i * wt_width);
    }
    
    //match bins for wt with jet pt included for e3c and Xbins for wt_res_e3c
    Double_t wt_from_e3c = -6;
    Double_t wt_to_e3c = 0;
    Int_t wt_bins_e3c = 400;
    Double_t wt_width_e3c = (wt_to_e3c-wt_from_e3c)/wt_bins_e3c;
    Double_t wt_new_bins_e3c[401] = {};
    for (Int_t i = 0; i <= wt_bins_e3c; i++)
    {
        wt_new_bins_e3c[i] = TMath::Power(10.0, wt_from_e3c + i * wt_width_e3c);
    }

    //for histograms to optimize unfolding bins 
    Double_t fromFineBins = -4;
    Double_t toFineBins = 0;
    Int_t binsFineBins = 50;
    Double_t widthFineBins = (toFineBins-fromFineBins)/binsFineBins;
    Double_t fine_bins[51] = {};
    for (Int_t i = 0; i <= binsFineBins; i++)
    {
        fine_bins[i] = TMath::Power(10.0, fromFineBins + i * widthFineBins);
    }

    Double_t from_constt = 40;
    Double_t to_constt = 200;
    Int_t bins_constt = 16;
    Double_t width_constt = (to_constt-from_constt)/bins_constt;
    Double_t tbins[17] = {};
    for (Int_t i = 0; i <= bins_constt; i++)
    {
        tbins[i] = (from_constt + i * width_constt);
    }

    Double_t from_const_reco = 40;
    Double_t to_const_reco = 100;
    Int_t bins_const_reco = 6;
    Double_t width_const_reco = (to_const_reco-from_const_reco)/bins_const_reco;
    Double_t xbins[9] = {};
    for (Int_t i = 0; i <= bins_const_reco; i++)
    {
        xbins[i] = (from_const_reco + i * width_const_reco);
    }
    xbins[7] = 120; 
    xbins[8] = 140; 

    double dRbins[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
        0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
    
    double dRbinst[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
        0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
    
    
    double wtbins[] = {
    0.00007, 0.0001, 0.000346737, 0.000912011, 0.00138038, 0.00165959, 0.00190546,
    0.00288403, 0.00346737, 0.00457088, 0.0060256, 0.00794328, 0.0104713, 
    0.0138038, 0.018197, 0.025, 0.0288403, 0.0346737, 0.0457088, 0.060256, 
    0.0794328, 0.104713, 0.138038, 0.18197, 0.25
};
    
    
    double wtbinst[] = {
    0.000001, 0.00005, 0.0001, 0.000346737, 0.0005, 0.000912011, 0.001, 0.00138038, 
    0.00165959, 0.00190546, 0.00288403, 0.00346737, 0.00457088, 0.0060256, 
    0.00794328, 0.00912011, 0.0104713, 0.0138038, 0.0165959, 0.018197, 
    0.0190546, 0.025, 0.0288403, 0.0346737, 0.0457088, 0.060256, 
    0.0794328, 0.104713, 0.138038, 0.18197, 0.25};

    double wtbins_e3c[] = {
        1e-06, 1.51356e-06, 5.62341e-06,1.01158e-05,1.12202e-05, 1.24451e-05,  1.38038e-05, 1.53109e-05,
        1.69824e-05,1.88365e-05, 
        2.0893e-05, 2.31739e-05, 2.5704e-05, 2.85102e-05, 3.16228e-05, 3.50752e-05,
        3.89045e-05, 4.31519e-05, 
        4.7863e-05, 5.30884e-05, 5.88844e-05, 6.53131e-05, 7.24436e-05, 8.03526e-05,
        8.91251e-05,  9.88553e-05, 0.000102329,0.000121619,0.000251189, 0.000501187,
        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158,  0.0108393,  0.0128825,  0.0153109,
        0.0175792, 0.0201837,
        0.0363078, 0.0462381, 0.0512861,0.0724436,
        0.0988553,
        0.1485};
    
    
    double wtbinst_e3c[] = {
        0.0001,0.0005,0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
        0.0288403,
        0.0346737,
        0.0457088, 0.060256,
        0.0794328,
        0.104713,0.138038, 0.18197,0.26
    };

    int nWtbins_e3c = sizeof(wtbins_e3c) / sizeof(wtbins_e3c[0]) - 1;
    int nWtbins = sizeof(wtbins) / sizeof(wtbins[0]) - 1;
    int ndRbins = sizeof(dRbins) / sizeof(dRbins[0]) - 1;
    int nJetPtbins = sizeof(xbins) / sizeof(xbins[0]) - 1;
    
    int nWtbinst_e3c = sizeof(wtbinst_e3c) / sizeof(wtbinst_e3c[0]) - 1;
    int nWtbinst = sizeof(wtbinst) / sizeof(wtbinst[0]) - 1;
    int ndRbinst = sizeof(dRbinst) / sizeof(dRbinst[0]) - 1;
    int nJetPtbinst = sizeof(tbins) / sizeof(tbins[0]) - 1;
    
    if(fCout){cout<<"HERE BINS ARE DECLARED!!!!!!!!!!!!!"<<endl;}
// /////////////////////////////// DATA HISTOGRAMS //////////////////////////////////////////
    if(!fDoEmbedding){
    jet_pt_hist = new TH1D("jet_pt_hist", "Jet Pt", 13,10,140);
    fOutput->Add(jet_pt_hist);

    pt_tru = new TH1D("jet_pt_tru_hist", "Jet Pt", 13,10,140);
    fOutput->Add(pt_tru);
    
    EEC_pt_hist = new TH2F("EEC_pt_hist", "EEC and jet_pt 2D", ndRbins, dRbins, 13,10,140);
    fOutput->Add(EEC_pt_hist);
    
    E3C_pt_hist = new TH2F("E3C_pt_hist", "EEEC and jet_pt 2D", ndRbins, dRbins, 13,10,140);
    fOutput->Add(E3C_pt_hist);

    delta_pt_cone = new TH1F("del_pt_cone","del_pt_cone",1000,-100,100); //all cones
    fOutput->Add(delta_pt_cone);

    delta_pt_coneBias1 = new TH1F("del_pt_cone_1","del_pt_cone_1Gev",1000,-100,100); //cones should at least have one 1 GeV particle
    fOutput->Add(delta_pt_coneBias1);

    delta_pt_coneBias2 = new TH1F("del_pt_cone_5Gev","del_pt_cone_5Gev",1000,-100,100); //cones should at least have one 5 GeV partilce
    fOutput->Add(delta_pt_coneBias2);

    delta_pt_coneBias3 = new TH1F("del_pt_cone_fMinPtConstMinus2","del_pt_cone_fMinPtConstMinus2",1000,-100,100); //cones should at least have one 5 GeV partilce
    fOutput->Add(delta_pt_coneBias3);

    delta_pt_coneBias4 = new TH1F("del_pt_cone_fMinPtConstPlus2","del_pt_cone_fMinPtConstPlus2",1000,-100,100); //cones should at least have one 5 GeV partilce
    fOutput->Add(delta_pt_coneBias4);

    delta_pt_ENCcone = new TH1F("del_pt_ENC_cone","del_pt_ENC_cone",1000,- 100,100); //cones should have all greater than 1 gev particles
    fOutput->Add(delta_pt_ENCcone);

    if(fCout) cout<<"###########------subtraction histograms for EEC data----########"<<endl;
  if(ifeec){
    h_MB1_dat = new TH2F("h_MB1_dat", "h_MB1_dat", 13, new_bins_const, ndRbins, dRbins);//min bias data
    fOutput->Add(h_MB1_dat);
    h_JMB_dat = new TH2F("h_JMB_dat", "h_JMB_dat", 13, new_bins_const, ndRbins, dRbins);//Jet*MINBIAS data
    fOutput->Add(h_JMB_dat);
    h_MB1MB2_dat = new TH2F("h_MB1MB2_dat", "h_MB1MB2_dat", 13, new_bins_const, ndRbins, dRbins);//Jet*MINBIAS data
    fOutput->Add(h_MB1MB2_dat);

    h3Jet_deltaR_MB1_dat = new TH3F("h3Jet_deltaR_MB1_dat", "h3Jet_deltaR_MB1_dat", nJetPtbins, xbins, ndRbins, dRbins,nWtbins,wtbins);//Jet*MINBIAS
    fOutput->Add(h3Jet_deltaR_MB1_dat);
    h3Jet_deltaR_JMB_dat = new TH3F("h3Jet_deltaR_JMB_dat", "h3Jet_deltaR_JMB_dat", nJetPtbins, xbins, ndRbins, dRbins,nWtbins,wtbins);//Jet*MINBIAS
    fOutput->Add(h3Jet_deltaR_JMB_dat);
    h3Jet_deltaR_MB1MB2_dat = new TH3F("h3Jet_deltaR_MB1MB2_dat", "h3_deltaR_MB1MB2_dat",nJetPtbins, xbins, ndRbins, dRbins,nWtbins,wtbins);//Jet*MINBIAS data
    fOutput->Add(h3Jet_deltaR_MB1MB2_dat);
  }
   if(ife3c){
   if(fCout) cout<<"###########------subtraction histograms for E3C data----########"<<endl;
    h_MB1MB1MB1_dat = new TH2F("h_MB1MB1MB1_dat", "h_MB1MB1MB1_dat",13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_MB1MB1MB1_dat);
    h_JJMB_dat = new TH2F("h_JJMB_dat", "h_JJMB_dat", 13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_JJMB_dat);
    h_JMB1MB2_dat = new TH2F("h_JMB1MB2_dat", "h_JMB1MB2_dat", 13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_JMB1MB2_dat);
    h_JMBMB_dat = new TH2F("h_JMBMB_dat", "h_JMBMB_dat",13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_JMBMB_dat);
    h_MB1MB2MB2_dat = new TH2F("h_MB1MB2MB2_dat", "h_MB1MB2MB2_dat",13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_MB1MB2MB2_dat);
    h_MB1MB1MB2_dat = new TH2F("h_MB1MB1MB2_dat", "h_MB1MB1MB2_dat",13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_MB1MB1MB2_dat);
    h_MB1MB2MB3_dat = new TH2F("h_MB1MB2MB3_dat", "h_MB1MB2MB3_dat",13, new_bins_const,ndRbins, dRbins);
    fOutput->Add(h_MB1MB2MB3_dat);

    h3_MB1MB1MB1_dat = new TH3F("h3_MB1MB1MB1_dat", "h3_MB1MB1MB1_dat",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB1_dat);
    h3_JJMB_dat = new TH3F("h3_JJMB_dat", "h3_JJMB_dat", nJetPtbins, xbins, ndRbins, dRbins,nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JJMB_dat);
    h3_JMB1MB2_dat = new TH3F("h3_JMB1MB2_dat", "h3_JMB1MB2_dat", nJetPtbins, xbins, ndRbins, dRbins,nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMB1MB2_dat);
    h3_JMBMB_dat = new TH3F("h3_JMBMB_dat", "h3_JMBMB_dat",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMBMB_dat);
    h3_MB1MB1MB2_dat = new TH3F("h3_MB1MB1MB2_dat", "h3_MB1MB1MB2_dat",nJetPtbins, xbins, ndRbins, dRbins,nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB2_dat);
    h3_MB1MB2MB2_dat = new TH3F("h3_MB1MB2MB2_dat", "h3_MB1MB2MB2_dat",nJetPtbins, xbins, ndRbins, dRbins,nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB2_dat);
    h3_MB1MB2MB3_dat = new TH3F("h3_MB1MB2MB3_dat", "h3_MB1MB2MB3_dat",nJetPtbins, xbins, ndRbins, dRbins,nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB3_dat);
   }

    }
    // //////////////////////////////////////////////////////////////////////////////////////////
   if(fDoEmbedding){
    hJet_num_um = new TH1F("hJet_num_um", "hJet_num_um",1,0,1);
    fOutput->Add(hJet_num_um);
    hJet_num_m = new TH1F("hJet_num_m", "hJet_num_m",1,0,1);
    fOutput->Add(hJet_num_m);
    h_jetpt_m = new TH1F("h_jetptM","h_jetptM",13,10,140);
    fOutput->Add(h_jetpt_m);
    h_jetpt_m_tru = new TH2F("h_jetptM_tru","h_jetptM_tru",13,10,140,13,10,140);
    fOutput->Add(h_jetpt_m_tru);
    h_jetpt_um = new TH1F("h_jetptUM","h_jetptUM",13,10,140);
    fOutput->Add(h_jetpt_um);
    h_dpt = new TH2F ("h_dpt","h_dpt",240,-60,60,13,10,140);
    fOutput->Add(h_dpt);
   }

   if(ifeec && fDoEmbedding)
   {
   if(fCout)cout<<"#####################---Declaring Embedding Histograms EEC----####################"<<endl;
    //matched subtracted jet
    h_MJ = new TH3F("h_MJ", "h_MJ", 13, new_bins_const,13, new_bins_const,21, new_bins);//all tracks fake + real
    fOutput->Add(h_MJ);
    h_MJ0 = new TH3F("h_MJ0", "h_MJ0",13, new_bins_const,13, new_bins_const,21, new_bins);//0 real tracks
    fOutput->Add(h_MJ0);
    h_MJ1 = new TH3F("h_MJ1", "h_MJ1", 13, new_bins_const,13, new_bins_const,21, new_bins);//1 real track
    fOutput->Add(h_MJ1);
    h_MJ2 = new TH3F("h_MJ2", "h_MJ2", 13, new_bins_const,13, new_bins_const,21, new_bins);//2 real tracks
    fOutput->Add(h_MJ2);

    h_MJ_m = new TH3F("h_MJ_m", "h_MJ_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//all tracks fake + real
    fOutput->Add(h_MJ_m);
    h_MJ0_m = new TH3F("h_MJ0_m", "h_MJ0_m",13, new_bins_const,13, new_bins_const,21, new_bins);//0 real tracks
    fOutput->Add(h_MJ0_m);
    h_MJ1_m = new TH3F("h_MJ1_m", "h_MJ1_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//1 real track
    fOutput->Add(h_MJ1_m);
    h_MJ2_m = new TH3F("h_MJ2_m", "h_MJ2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//2 real tracks
    fOutput->Add(h_MJ2_m);
    
    h_MJ_tru_m = new TH3F("h_MJ_tru_m", "h_MJ_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//all tracks fake + real
    fOutput->Add(h_MJ_tru_m);
    h_MJ0_tru_m = new TH3F("h_MJ0_tru_m", "h_MJ0_tru_m",13, new_bins_const,13, new_bins_const,21, new_bins);//0 real tracks
    fOutput->Add(h_MJ0_tru_m);
    h_MJ1_tru_m = new TH3F("h_MJ1_tru_m", "h_MJ1_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//1 real track
    fOutput->Add(h_MJ1_tru_m);
    h_MJ2_tru_m = new TH3F("h_MJ2_tru_m", "h_MJ2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//2 real tracks
    fOutput->Add(h_MJ2_tru_m);

    h_MJ_um = new TH3F("h_MJ_um", "h_MJ_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//all tracks fake + real
    fOutput->Add(h_MJ_um);
    h_MJ0_um = new TH3F("h_MJ0_um", "h_MJ0_um",13, new_bins_const,13, new_bins_const,21, new_bins);//0 real tracks
    fOutput->Add(h_MJ0_um);
    h_MJ1_um = new TH3F("h_MJ1_um", "h_MJ1_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//1 real track
    fOutput->Add(h_MJ1_um);
    h_MJ2_um = new TH3F("h_MJ2_um", "h_MJ2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//2 real tracks
    fOutput->Add(h_MJ2_um);
    
  //for MB1
    h_MB1 = new TH3F("h_MB1", "h_MB1", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1);
    h_MB1_m = new TH3F("h_MB1_m", "h_MB1_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1_m);
    h_MB1_tru_m = new TH3F("h_MB1_tru_m", "h_MB1_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1_tru_m);
    h_MB1_um = new TH3F("h_MB1_um", "h_MB1_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1_um);

  //for JMB
    h_JMB = new TH3F("h_JMB", "h_JMB", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_JMB);
    h_JMB_m = new TH3F("h_JMB_m", "h_JMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_JMB_m);
    h_JMB_tru_m = new TH3F("h_JMB_tru_m", "h_JMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_JMB_tru_m);
    h_JMB_um = new TH3F("h_JMB_um", "h_JMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_JMB_um);
   
    h_SMB = new TH3F("h_SMB", "h_SMB", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_SMB);
    h_SMB_m = new TH3F("h_SMB_m", "h_SMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_SMB_m);
    h_SMB_tru_m = new TH3F("h_SMB_tru_m", "h_SMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_SMB_tru_m);
    h_SMB_um = new TH3F("h_SMB_um", "h_SMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_SMB_um);

    h_BMB = new TH3F("h_BMB", "h_BMB", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_BMB);
    h_BMB_m = new TH3F("h_BMB_m", "h_BMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_BMB_m);
    h_BMB_tru_m = new TH3F("h_BMB_tru_m", "h_BMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_BMB_tru_m);
    h_BMB_um = new TH3F("h_BMB_um", "h_BMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_BMB_um);

  //for MB1MB2
    h_MB1MB2 = new TH3F("h_MB1MB2", "h_MB1MB2", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1MB2);
    h_MB1MB2_m = new TH3F("h_MB1MB2_m", "h_MB1MB2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1MB2_m);
    h_MB1MB2_tru_m = new TH3F("h_MB1MB2_tru_m", "h_MB1MB2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1MB2_tru_m);
    h_MB1MB2_um = new TH3F("h_MB1MB2_um", "h_MB1MB2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);//min bias
    fOutput->Add(h_MB1MB2_um);

//////////////////////// HISTOGRAMS FOR 3D Subtraction and Unfolding///////////////////////////

if(fCout)cout<<"#####################---Declaring 3D Embedding Histograms EEC----####################"<<endl;
   //For same jet 
    h3Jet_deltaR_MJ = new TH3F("h3Jet_deltaR_MJ", "h3Jet_deltaR_MJ", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//all tracks fake + real
    fOutput->Add(h3Jet_deltaR_MJ);
    h3Jet_deltaR_MJ0 = new TH3F("h3Jet_deltaR_MJ0", "h3Jet_deltaR_MJ0",nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//0 real tracks
    fOutput->Add(h3Jet_deltaR_MJ0);
    h3Jet_deltaR_MJ1 = new TH3F("h3Jet_deltaR_MJ1", "h3Jet_deltaR_MJ1", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//1 real track
    fOutput->Add(h3Jet_deltaR_MJ1);
    h3Jet_deltaR_MJ2 = new TH3F("h3Jet_deltaR_MJ2", "h3Jet_deltaR_MJ2", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//2 real tracks
    fOutput->Add(h3Jet_deltaR_MJ2);

    h3Jet_deltaR_MJ_m = new TH3F("h3Jet_deltaR_MJ_m", "h3Jet_deltaR_MJ_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//all tracks fake + real
    fOutput->Add(h3Jet_deltaR_MJ_m);
    h3Jet_deltaR_MJ0_m = new TH3F("h3Jet_deltaR_MJ0_m", "h3Jet_deltaR_MJ0_m",nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//0 real tracks
    fOutput->Add(h3Jet_deltaR_MJ0_m);
    h3Jet_deltaR_MJ1_m = new TH3F("h3Jet_deltaR_MJ1_m", "h3Jet_deltaR_MJ1_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//1 real track
    fOutput->Add(h3Jet_deltaR_MJ1_m);
    h3Jet_deltaR_MJ2_m = new TH3F("h3Jet_deltaR_MJ2_m", "h3Jet_deltaR_MJ2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//2 real tracks
    fOutput->Add(h3Jet_deltaR_MJ2_m);

    h3Jet_deltaR_MJ_tru_m = new TH3F("h3Jet_deltaR_MJ_tru_m", "h3Jet_deltaR_MJ_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//all tracks fake + real
    fOutput->Add(h3Jet_deltaR_MJ_tru_m);
    h3Jet_deltaR_MJ0_tru_m = new TH3F("h3Jet_deltaR_MJ0_tru_m", "h3Jet_deltaR_MJ0_tru_m",nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//0 real tracks
    fOutput->Add(h3Jet_deltaR_MJ0_tru_m);
    h3Jet_deltaR_MJ1_tru_m = new TH3F("h3Jet_deltaR_MJ1_tru_m", "h3Jet_deltaR_MJ1_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//1 real track
    fOutput->Add(h3Jet_deltaR_MJ1_tru_m);
    h3Jet_deltaR_MJ2_tru_m = new TH3F("h3Jet_deltaR_MJ2_tru_m", "h3Jet_deltaR_MJ2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//2 real tracks
    fOutput->Add(h3Jet_deltaR_MJ2_tru_m);

    h3Jet_deltaR_MJ_um = new TH3F("h3Jet_deltaR_MJ_um", "h3Jet_deltaR_MJ_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//all tracks fake + real
    fOutput->Add(h3Jet_deltaR_MJ_um);
    h3Jet_deltaR_MJ0_um = new TH3F("h3Jet_deltaR_MJ0_um", "h3Jet_deltaR_MJ0_um",nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//0 real tracks
    fOutput->Add(h3Jet_deltaR_MJ0_um);
    h3Jet_deltaR_MJ1_um = new TH3F("h3Jet_deltaR_MJ1_um", "h3Jet_deltaR_MJ1_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//1 real track
    fOutput->Add(h3Jet_deltaR_MJ1_um);
    h3Jet_deltaR_MJ2_um = new TH3F("h3Jet_deltaR_MJ2_um", "h3Jet_deltaR_MJ2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//2 real tracks
    fOutput->Add(h3Jet_deltaR_MJ2_um);

//for MB1
    h3Jet_deltaR_MB1 = new TH3F("h3Jet_deltaR_MB1", "h3Jet_deltaR_MB1", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins); 
    fOutput->Add(h3Jet_deltaR_MB1);
    h3Jet_deltaR_MB1_m = new TH3F("h3Jet_deltaR_MB1_m", "h3Jet_deltaR_MB1_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);  
    fOutput->Add(h3Jet_deltaR_MB1_m);
    h3Jet_deltaR_MB1_tru_m = new TH3F("h3Jet_deltaR_MB1_tru_m", "h3Jet_deltaR_MB1_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);   
    fOutput->Add(h3Jet_deltaR_MB1_tru_m);
    h3Jet_deltaR_MB1_um = new TH3F("h3Jet_deltaR_MB1_um", "h3Jet_deltaR_MB1_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);  
    fOutput->Add(h3Jet_deltaR_MB1_um); 
//for JMB
    h3Jet_deltaR_JMB = new TH3F("h3Jet_deltaR_JMB", "h3Jet_deltaR_JMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//Jet*MINBIAS
    fOutput->Add(h3Jet_deltaR_JMB);
    h3Jet_deltaR_JMB_m = new TH3F("h3Jet_deltaR_JMB_m", "h3Jet_deltaR_JMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//Jet*MINBIAS
    fOutput->Add(h3Jet_deltaR_JMB_m);
    h3Jet_deltaR_JMB_tru_m = new TH3F("h3Jet_deltaR_JMB_tru_m", "h3Jet_deltaR_JMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);
    fOutput->Add(h3Jet_deltaR_JMB_tru_m);
    h3Jet_deltaR_JMB_um = new TH3F("h3Jet_deltaR_JMB_um", "h3Jet_deltaR_JMB_tru_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);
    fOutput->Add(h3Jet_deltaR_JMB_um);
//for SMB
    h3Jet_deltaR_SMB = new TH3F("h3Jet_deltaR_SMB", "h3Jet_deltaR_SMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_SMB);
    h3Jet_deltaR_SMB_m = new TH3F("h3Jet_deltaR_SMB_m", "h3Jet_deltaR_SMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_SMB_m);
    h3Jet_deltaR_SMB_tru_m = new TH3F("h3Jet_deltaR_SMB_tru_m", "h3Jet_deltaR_SMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_SMB_tru_m);
    h3Jet_deltaR_SMB_um = new TH3F("h3Jet_deltaR_SMB_um", "h3Jet_deltaR_SMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_SMB_um);
//for BMB 
    h3Jet_deltaR_BMB = new TH3F("h3Jet_deltaR_BMB", "h3Jet_deltaR_BMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_BMB);
    h3Jet_deltaR_BMB_m = new TH3F("h3Jet_deltaR_BMB_m", "h3Jet_deltaR_BMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_BMB_m);
    h3Jet_deltaR_BMB_tru_m = new TH3F("h3Jet_deltaR_BMB_tru_m", "h3Jet_deltaR_BMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_BMB_tru_m);
    h3Jet_deltaR_BMB_um = new TH3F("h3Jet_deltaR_BMB_um", "h3Jet_deltaR_BMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_BMB_um);
//for MB1MB2
    h3Jet_deltaR_MB1MB2 = new TH3F("h3Jet_deltaR_MB1MB2", "h3Jet_deltaR_MB1MB2", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_MB1MB2);
    h3Jet_deltaR_MB1MB2_m = new TH3F("h3Jet_deltaR_MB1MB2_m", "h3Jet_deltaR_MB1MB2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_MB1MB2_m);
    h3Jet_deltaR_MB1MB2_tru_m = new TH3F("h3Jet_deltaR_MB1MB2_tru_m", "h3Jet_deltaR_MB1MB2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_MB1MB2_tru_m);
    h3Jet_deltaR_MB1MB2_um = new TH3F("h3Jet_deltaR_MB1MB2_um", "h3Jet_deltaR_MB1MB2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);     
    fOutput->Add(h3Jet_deltaR_MB1MB2_um);


    if(ifMinPtHist){
    if(fCout){cout<<"################!!!!!!!!!!!!Creating histograms with truth pT min cut for embedded jets##############!!!!!!!!!!"<<endl;}
    h3Jet_deltaR_MJ_minpt = new TH3F("h3Jet_deltaR_MJ_minpt", "h3Jet_deltaR_MJ_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//all tracks fake + real
    fOutput->Add(h3Jet_deltaR_MJ_minpt);
    h3Jet_deltaR_MJ0_minpt = new TH3F("h3Jet_deltaR_MJ0_minpt", "h3Jet_deltaR_MJ0_minpt",nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//0 real tracks
    fOutput->Add(h3Jet_deltaR_MJ0_minpt);
    h3Jet_deltaR_MJ1_minpt = new TH3F("h3Jet_deltaR_MJ1_minpt", "h3Jet_deltaR_MJ1_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//1 real track
    fOutput->Add(h3Jet_deltaR_MJ1_minpt);
    h3Jet_deltaR_MJ2_minpt = new TH3F("h3Jet_deltaR_MJ2_minpt", "h3Jet_deltaR_MJ2_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//2 real tracks
    fOutput->Add(h3Jet_deltaR_MJ2_minpt);
    h3Jet_deltaR_MB1_minpt = new TH3F("h3Jet_deltaR_MB1_minpt", "h3Jet_deltaR_MB1_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//min bias    //Min bias 1
    fOutput->Add(h3Jet_deltaR_MB1_minpt);
    h3Jet_deltaR_JMB_minpt = new TH3F("h3Jet_deltaR_JMB_minpt", "h3Jet_deltaR_JMB_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//Jet*MINBIAS
    fOutput->Add(h3Jet_deltaR_JMB_minpt);
    h3Jet_deltaR_SMB_minpt = new TH3F("h3Jet_deltaR_SMB_minpt", "h3Jet_deltaR_SMB_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//JetSignal*MINBIAS
    fOutput->Add(h3Jet_deltaR_SMB_minpt);
    h3Jet_deltaR_BMB_minpt = new TH3F("h3Jet_deltaR_BMB_minpt", "h3Jet_deltaR_BMB_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//JetBkg*MINBIAS
    fOutput->Add(h3Jet_deltaR_BMB_minpt);
    h3Jet_deltaR_MB1MB2_minpt = new TH3F("h3Jet_deltaR_MB1MB2_minpt", "h3Jet_deltaR_MB1MB2_minpt", nJetPtbins, xbins, ndRbins, dRbins, nWtbins,wtbins);//MinBias1*MinBias2
    fOutput->Add(h3Jet_deltaR_MB1MB2_minpt);
    h3Jet_deltaR_MJ_tru_minpt = new TH3F("h3Jet_deltaR_MJ_tru_minpt", "h3Jet_deltaR_MJ_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//all tracks fake + real
    fOutput->Add(h3Jet_deltaR_MJ_tru_minpt);
    h3Jet_deltaR_MJ0_tru_minpt = new TH3F("h3Jet_deltaR_MJ0_tru_minpt", "h3Jet_deltaR_MJ0_tru_minpt",nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//0 real tracks
    fOutput->Add(h3Jet_deltaR_MJ0_tru_minpt);
    h3Jet_deltaR_MJ1_tru_minpt = new TH3F("h3Jet_deltaR_MJ1_tru_minpt", "h3Jet_deltaR_MJ1_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//1 real track
    fOutput->Add(h3Jet_deltaR_MJ1_tru_minpt);
    h3Jet_deltaR_MJ2_tru_minpt = new TH3F("h3Jet_deltaR_MJ2_tru_minpt", "h3Jet_deltaR_MJ2_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//2 real tracks
    fOutput->Add(h3Jet_deltaR_MJ2_tru_minpt);
    h3Jet_deltaR_MB1_tru_minpt = new TH3F("h3Jet_deltaR_MB1_tru_minpt", "h3Jet_deltaR_MB1_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//min bias
    fOutput->Add(h3Jet_deltaR_MB1_tru_minpt);
    h3Jet_deltaR_JMB_tru_minpt = new TH3F("h3Jet_deltaR_JMB_tru_minpt", "h3Jet_deltaR_JMB_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//Jet*MINBIAS
    fOutput->Add(h3Jet_deltaR_JMB_tru_minpt);
    h3Jet_deltaR_SMB_tru_minpt = new TH3F("h3Jet_deltaR_SMB_tru_minpt", "h3Jet_deltaR_SMB_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//JetSignal*MINBIAS
    fOutput->Add(h3Jet_deltaR_SMB_tru_minpt);
    h3Jet_deltaR_BMB_tru_minpt = new TH3F("h3Jet_deltaR_BMB_tru_minpt", "h3Jet_deltaR_BMB_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//JetBkg*MINBIAS
    fOutput->Add(h3Jet_deltaR_BMB_tru_minpt);
    h3Jet_deltaR_MB1MB2_tru_minpt = new TH3F("h3Jet_deltaR_MB1MB2_tru_minpt", "h3Jet_deltaR_MB1MB2_tru_minpt", nJetPtbinst, tbins, ndRbinst, dRbinst, nWtbinst,wtbinst);//MinBias1*MinBias2
    fOutput->Add(h3Jet_deltaR_MB1MB2_tru_minpt);
    }
  }

    if(fCout){cout<<"#######################!!!!!!!!!!!!!!!All declarations for EEC successful!!!!!!!!!!!!!####################"<<endl;}
    if(fCout) cout<<"subtraction histograms for E3C"<<endl;
//     /////FOR E3C MIN BIAS SUBTRCTION////////
  if(ife3c && fDoEmbedding)
  {
     h_MJ_e3c = new TH3F("h_MJ_e3c", "h_MJ_e3c", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ_e3c);
     h_MJ0_e3c = new TH3F("h_MJ0_e3c", "h_MJ0_e3c", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ0_e3c);
     h_MJ1_e3c = new TH3F("h_MJ1_e3c", "h_MJ1_e3c",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ1_e3c);
     h_MJ2_e3c = new TH3F("h_MJ2_e3c", "h_MJ2_e3c",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ2_e3c);
     h_MJ3_e3c = new TH3F("h_MJ3_e3c", "h_MJ3_e3c", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ3_e3c);
     h_MJ_e3c_tru = new TH3F("h_MJ_e3c_tru", "h_MJ_e3c_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ_e3c_tru);
     h_MJ0_e3c_tru = new TH3F("h_MJ0_e3c_tru", "h_MJ0_e3c_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ0_e3c_tru);
     h_MJ1_e3c_tru = new TH3F("h_MJ1_e3c_tru", "h_MJ1_e3c_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ1_e3c_tru);
     h_MJ2_e3c_tru = new TH3F("h_MJ2_e3c_tru", "h_MJ2_e3c_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ2_e3c_tru);
     h_MJ3_e3c_tru = new TH3F("h_MJ3_e3c_tru", "h_MJ3_e3c_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MJ3_e3c_tru);

    if(fCout) cout<<"subtraction histograms for E3C 1"<<endl; 

     h_MB1MB1MB1 = new TH3F("h_MB1MB1MB1", "h_MB1MB1MB1",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB1MB1);
     h_MB1MB1MB1_tru = new TH3F("h_MB1MB1MB1_tru", "h_MB1MB1MB1_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB1MB1_tru);

    if(fCout) cout<<"subtraction histograms for E3C 2"<<endl; 

     h_JJMB = new TH3F("h_JJMB", "h_JJMB", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_JJMB);
     h_JJMB_tru = new TH3F("h_JJMB_tru", "h_JJMB_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_JJMB_tru);

     h_JMBMB = new TH3F("h_JMBMB", "h_JMBMB",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_JMBMB);
     h_JMBMB_tru = new TH3F("h_JMBMB_tru", "h_JMBMB_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_JMBMB_tru);

     h_MB1MB1MB2 = new TH3F("h_MB1MB1MB2", "h_MB1MB1MB2",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB1MB2);
     h_MB1MB1MB2_tru = new TH3F("h_MB1MB1MB2_tru", "h_MB1MB1MB2_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB1MB2_tru);

     h_MB1MB2MB2 = new TH3F("h_MB1MB2MB2", "h_MB1MB2MB2",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB2MB2);
     h_MB1MB2MB2_tru = new TH3F("h_MB1MB2MB2_tru", "h_MB1MB2MB2_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB2MB2_tru);

    if(fCout) cout<<"subtraction histograms for E3C 3"<<endl; 

     h_JMB1MB2 = new TH3F("h_JMB1MB2", "h_JMB1MB2", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_JMB1MB2);
     h_JMB1MB2_tru = new TH3F("h_JMB1MB2_tru", "h_JMB1MB2_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_JMB1MB2_tru);

     h_MB1MB2MB3 = new TH3F("h_MB1MB2MB3", "h_MB1MB2MB3",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB2MB3);
     h_MB1MB2MB3_tru = new TH3F("h_MB1MB2MB3_tru", "h_MB1MB2MB3_tru",13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_MB1MB2MB3_tru);

    if(fCout) cout<<"subtraction histograms for E3C 4"<<endl; 

     h_BMBMB = new TH3F("h_BMBMB", "h_BMBMB", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_BMBMB);
     h_SMBMB = new TH3F("h_SMBMB", "h_SMBMB", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_SMBMB);
     h_BMB1MB2 = new TH3F("h_BMB1MB2", "h_BMB1MB2", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_BMB1MB2);
     h_SMB1MB2 = new TH3F("h_SMB1MB2", "h_SMB1MB2", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_SMB1MB2);
     h_BBMB = new TH3F("h_BBMB", "h_BBMB", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_BBMB);
     h_SBMB = new TH3F("h_SBMB", "h_SBMB", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_SBMB);
     h_SSMB = new TH3F("h_SSMB", "h_SSMB", 13, new_bins_const,13, new_bins_const,21, new_bins);
     fOutput->Add(h_SSMB);

    h_BMBMB_tru = new TH3F("h_BMBMB_tru", "h_BMBMB_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMBMB_tru);
    h_SMBMB_tru = new TH3F("h_SMBMB_tru", "h_SMBMB_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMBMB_tru);
    h_BMB1MB2_tru = new TH3F("h_BMB1MB2_tru", "h_BMB1MB2_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMB1MB2_tru);
    h_BBMB_tru = new TH3F("h_BBMB_tru", "h_BBMB_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BBMB_tru);
    h_SBMB_tru = new TH3F("h_SBMB_tru", "h_SBMB_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SBMB_tru);
    h_SSMB_tru = new TH3F("h_SSMB_tru", "h_SSMB_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SSMB_tru);
    h_SMB1MB2_tru = new TH3F("h_SMB1MB2_tru", "h_SMB1MB2_tru", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMB1MB2_tru);

    if(fCout) cout<<"subtraction histograms for E3C matched jets "<<endl;
    hJet_deltaR_MJ_e3c_m = new TH3F("hJet_deltaR_MJ_e3c_m", "hJet_deltaR_MJ_e3c_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ_e3c_m);
    hJet_deltaR_MJ0_e3c_m = new TH3F("hJet_deltaR_MJ0_e3c_m", "hJet_deltaR_MJ0_e3c_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ0_e3c_m);
    hJet_deltaR_MJ1_e3c_m = new TH3F("hJet_deltaR_MJ1_e3c_m", "hJet_deltaR_MJ1_e3c_m",13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ1_e3c_m);
    hJet_deltaR_MJ2_e3c_m = new TH3F("hJet_deltaR_MJ2_e3c_m", "hJet_deltaR_MJ2_e3c_m",13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ2_e3c_m);
    hJet_deltaR_MJ3_e3c_m = new TH3F("hJet_deltaR_MJ3_e3c_m", "hJet_deltaR_MJ3_e3c_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ3_e3c_m);
    h_MB1MB1MB1_m = new TH3F("h_MB1MB1MB1_m", "h_MB1MB1MB1_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB1MB1_m);
    h_JJMB_m = new TH3F("h_JJMB_m", "h_JJMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JJMB_m);
    h_JMBMB_m = new TH3F("h_JMBMB_m", "h_JMBMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JMBMB_m);
    h_MB1MB1MB2_m = new TH3F("h_MB1MB1MB2_m", "h_MB1MB1MB2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB1MB2_m);
    h_MB1MB2MB2_m = new TH3F("h_MB1MB2MB2_m", "h_MB1MB2MB2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB2MB2_m);
    h_JMB1MB2_m = new TH3F("h_JMB1MB2_m", "h_JMB1MB2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JMB1MB2_m);
    h_MB1MB2MB3_m = new TH3F("h_MB1MB2MB3_m", "h_MB1MB2MB3_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB2MB3_m);
    h_BMBMB_m = new TH3F("h_BMBMB_m", "h_BMBMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMBMB_m);
    h_SMBMB_m = new TH3F("h_SMBMB_m", "h_SMBMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMBMB_m);
    h_BMB1MB2_m = new TH3F("h_BMB1MB2_m", "h_BMB1MB2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMB1MB2_m);
    h_SMB1MB2_m = new TH3F("h_SMB1MB2_m", "h_SMB1MB2_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMB1MB2_m);
    h_BBMB_m = new TH3F("h_BBMB_m", "h_BBMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BBMB_m);
    h_SBMB_m = new TH3F("h_SBMB_m", "h_SBMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SBMB_m);
    h_SSMB_m = new TH3F("h_SSMB_m", "h_SSMB_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SSMB_m);
    
    hJet_deltaR_MJ_e3c_tru_m = new TH3F("hJet_deltaR_MJ_e3c_tru_m", "hJet_deltaR_MJ_e3c_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ_e3c_tru_m);
    hJet_deltaR_MJ0_e3c_tru_m = new TH3F("hJet_deltaR_MJ0_e3c_tru_m", "hJet_deltaR_MJ0_e3c_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ0_e3c_tru_m);
    hJet_deltaR_MJ1_e3c_tru_m = new TH3F("hJet_deltaR_MJ1_e3c_tru_m", "hJet_deltaR_MJ1_e3c_tru_m",13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ1_e3c_tru_m);
    hJet_deltaR_MJ2_e3c_tru_m = new TH3F("hJet_deltaR_MJ2_e3c_tru_m", "hJet_deltaR_MJ2_e3c_tru_m",13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ2_e3c_tru_m);
    hJet_deltaR_MJ3_e3c_tru_m = new TH3F("hJet_deltaR_MJ3_e3c_tru_m", "hJet_deltaR_MJ3_e3c_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ3_e3c_tru_m);
    h_JJMB_tru_m = new TH3F("h_JJMB_tru_m", "h_JJMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JJMB_tru_m);
    h_JMBMB_tru_m = new TH3F("h_JMBMB_tru_m", "h_JMBMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JMBMB_tru_m);
    h_MB1MB1MB2_tru_m = new TH3F("h_MB1MB1MB2_tru_m", "h_MB1MB1MB2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB1MB2_tru_m);
    h_MB1MB2MB2_tru_m = new TH3F("h_MB1MB2MB2_tru_m", "h_MB1MB2MB2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB2MB2_tru_m);
    h_JMB1MB2_tru_m = new TH3F("h_JMB1MB2_tru_m", "h_JMB1MB2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JMB1MB2_tru_m);
    h_MB1MB2MB3_tru_m = new TH3F("h_MB1MB2MB3_tru_m", "h_MB1MB2MB3_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB2MB3_tru_m);
    h_BMBMB_tru_m = new TH3F("h_BMBMB_tru_m", "h_BMBMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMBMB_tru_m);
    h_SMBMB_tru_m = new TH3F("h_SMBMB_tru_m", "h_SMBMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMBMB_tru_m);
    h_BMB1MB2_tru_m = new TH3F("h_BMB1MB2_tru_m", "h_BMB1MB2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMB1MB2_tru_m);
    h_BBMB_tru_m = new TH3F("h_BBMB_tru_m", "h_BBMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BBMB_tru_m);
    h_SBMB_tru_m = new TH3F("h_SBMB_tru_m", "h_SBMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SBMB_tru_m);
    h_SSMB_tru_m = new TH3F("h_SSMB_tru_m", "h_SSMB_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SSMB_tru_m);
    h_SMB1MB2_tru_m = new TH3F("h_SMB1MB2_tru_m", "h_SMB1MB2_tru_m", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMB1MB2_tru_m);

  if(fCout) cout<<"subtraction histograms for E3C unmatched jets "<<endl;
    hJet_deltaR_MJ_e3c_um = new TH3F("hJet_deltaR_MJ_e3c_um", "hJet_deltaR_MJ_e3c_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ_e3c_um);
    hJet_deltaR_MJ0_e3c_um = new TH3F("hJet_deltaR_MJ0_e3c_um", "hJet_deltaR_MJ0_e3c_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ0_e3c_um);
    hJet_deltaR_MJ1_e3c_um = new TH3F("hJet_deltaR_MJ1_e3c_um", "hJet_deltaR_MJ1_e3c_um",13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ1_e3c_um);
    hJet_deltaR_MJ2_e3c_um = new TH3F("hJet_deltaR_MJ2_e3c_um", "hJet_deltaR_MJ2_e3c_um",13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ2_e3c_um);
    hJet_deltaR_MJ3_e3c_um = new TH3F("hJet_deltaR_MJ3_e3c_um", "hJet_deltaR_MJ3_e3c_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(hJet_deltaR_MJ3_e3c_um);
    h_MB1MB1MB1_um = new TH3F("h_MB1MB1MB1_um", "h_MB1MB1MB1_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB1MB1_um);
    h_JJMB_um = new TH3F("h_JJMB_um", "h_JJMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JJMB_um);
    h_JMBMB_um = new TH3F("h_JMBMB_um", "h_JMBMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JMBMB_um);
    h_MB1MB1MB2_um = new TH3F("h_MB1MB1MB2_um", "h_MB1MB1MB2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB1MB2_um);
    h_MB1MB2MB2_um = new TH3F("h_MB1MB2MB2_um", "h_MB1MB2MB2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB2MB2_um);
    h_JMB1MB2_um = new TH3F("h_JMB1MB2_um", "h_JMB1MB2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_JMB1MB2_um);
    h_MB1MB2MB3_um = new TH3F("h_MB1MB2MB3_um", "h_MB1MB2MB3_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_MB1MB2MB3_um);
    h_BMBMB_um = new TH3F("h_BMBMB_um", "h_BMBMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMBMB_um);
    h_SMBMB_um = new TH3F("h_SMBMB_um", "h_SMBMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMBMB_um);
    h_BMB1MB2_um = new TH3F("h_BMB1MB2_um", "h_BMB1MB2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BMB1MB2_um);
    h_SMB1MB2_um = new TH3F("h_SMB1MB2_um", "h_SMB1MB2_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SMB1MB2_um);
    h_BBMB_um = new TH3F("h_BBMB_um", "h_BBMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_BBMB_um);
    h_SBMB_um = new TH3F("h_SBMB_um", "h_SBMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SBMB_um);
    h_SSMB_um = new TH3F("h_SSMB_um", "h_SSMB_um", 13, new_bins_const,13, new_bins_const,21, new_bins);
    fOutput->Add(h_SSMB_um);

  if(fCout) cout<<"subtraction histograms for E3C unfolding bins 3D"<<endl;
    h3Jet_deltaR_MJ_e3c = new TH3F("h3Jet_deltaR_MJ_e3c", "h3Jet_deltaR_MJ_e3c", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ_e3c);
    h3Jet_deltaR_MJ0_e3c = new TH3F("h3Jet_deltaR_MJ0_e3c", "h3Jet_deltaR_MJ0_e3c", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ0_e3c);
    h3Jet_deltaR_MJ1_e3c = new TH3F("h3Jet_deltaR_MJ1_e3c", "h3Jet_deltaR_MJ1_e3c",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ1_e3c);
    h3Jet_deltaR_MJ2_e3c = new TH3F("h3Jet_deltaR_MJ2_e3c", "h3Jet_deltaR_MJ2_e3c",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ2_e3c);
    h3Jet_deltaR_MJ3_e3c = new TH3F("h3Jet_deltaR_MJ3_e3c", "h3Jet_deltaR_MJ3_e3c", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ3_e3c);
    
    h3_MB1MB1MB1 = new TH3F("h3_MB1MB1MB1", "h3_MB1MB1MB1",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB1);
    h3_JJMB = new TH3F("h3_JJMB", "h3_JJMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JJMB);
    h3_JMBMB = new TH3F("h3_JMBMB", "h3_JMBMB",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMBMB);
    h3_MB1MB1MB2 = new TH3F("h3_MB1MB1MB2", "h3_MB1MB1MB2",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB2);
    h3_MB1MB2MB2 = new TH3F("h3_MB1MB2MB2", "h3_MB1MB2MB2",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB2);
    h3_JMB1MB2 = new TH3F("h3_JMB1MB2", "h3_JMB1MB2", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMB1MB2);
    h3_MB1MB2MB3 = new TH3F("h3_MB1MB2MB3", "h3_MB1MB2MB3",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB3);
    h3_BMBMB = new TH3F("h3_BMBMB", "h3_BMBMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BMBMB);
    h3_SMBMB = new TH3F("h3_SMBMB", "h3_SMBMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SMBMB);
    h3_BMB1MB2 = new TH3F("h3_BMB1MB2", "h3_BMB1MB2", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BMB1MB2);
    h3_SMB1MB2 = new TH3F("h3_SMB1MB2", "h3_SMB1MB2", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SMB1MB2);
    h3_BBMB = new TH3F("h3_BBMB", "h3_BBMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BBMB);
    h3_SBMB = new TH3F("h3_SBMB", "h3_SBMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SBMB);
    h3_SSMB = new TH3F("h3_SSMB", "h3_SSMB", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SSMB);

  if(fCout) cout<<"matched jets only"<<endl;
    h3Jet_deltaR_MJ_e3c_m = new TH3F("h3Jet_deltaR_MJ_e3c_m", "h3Jet_deltaR_MJ_e3c_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ_e3c_m);
    h3Jet_deltaR_MJ0_e3c_m = new TH3F("h3Jet_deltaR_MJ0_e3c_m", "h3Jet_deltaR_MJ0_e3c_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ0_e3c_m);
    h3Jet_deltaR_MJ1_e3c_m = new TH3F("h3Jet_deltaR_MJ1_e3c_m", "h3Jet_deltaR_MJ1_e3c_m",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ1_e3c_m);
    h3Jet_deltaR_MJ2_e3c_m = new TH3F("h3Jet_deltaR_MJ2_e3c_m", "h3Jet_deltaR_MJ2_e3c_m",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ2_e3c_m);
    h3Jet_deltaR_MJ3_e3c_m = new TH3F("h3Jet_deltaR_MJ3_e3c_m", "h3Jet_deltaR_MJ3_e3c_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ3_e3c_m);
    h3_MB1MB1MB1_m = new TH3F("h3_MB1MB1MB1_m", "h3_MB1MB1MB1_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB1_m);
    h3_JJMB_m = new TH3F("h3_JJMB_m", "h3_JJMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JJMB_m);
    h3_JMBMB_m = new TH3F("h3_JMBMB_m", "h3_JMBMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMBMB_m);
    h3_MB1MB1MB2_m = new TH3F("h3_MB1MB1MB2_m", "h3_MB1MB1MB2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB2_m);
    h3_MB1MB2MB2_m = new TH3F("h3_MB1MB2MB2_m", "h3_MB1MB2MB2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB2_m);
    h3_JMB1MB2_m = new TH3F("h3_JMB1MB2_m", "h3_JMB1MB2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMB1MB2_m);
    h3_MB1MB2MB3_m = new TH3F("h3_MB1MB2MB3_m", "h3_MB1MB2MB3_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB3_m);
    h3_BMBMB_m = new TH3F("h3_BMBMB_m", "h3_BMBMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BMBMB_m);
    h3_SMBMB_m = new TH3F("h3_SMBMB_m", "h3_SMBMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SMBMB_m);
    h3_BMB1MB2_m = new TH3F("h3_BMB1MB2_m", "h3_BMB1MB2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BMB1MB2_m);
    h3_SMB1MB2_m = new TH3F("h3_SMB1MB2_m", "h3_SMB1MB2_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SMB1MB2_m);
    h3_BBMB_m = new TH3F("h3_BBMB_m", "h3_BBMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BBMB_m);
    h3_SBMB_m = new TH3F("h3_SBMB_m", "h3_SBMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SBMB_m);
    h3_SSMB_m = new TH3F("h3_SSMB_m", "h3_SSMB_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SSMB_m);
    
    h3Jet_deltaR_MJ_e3c_tru_m = new TH3F("h3Jet_deltaR_MJ_e3c_tru_m", "h3Jet_deltaR_MJ_e3c_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3Jet_deltaR_MJ_e3c_tru_m);
    h3Jet_deltaR_MJ0_e3c_tru_m = new TH3F("h3Jet_deltaR_MJ0_e3c_tru_m", "h3Jet_deltaR_MJ0_e3c_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3Jet_deltaR_MJ0_e3c_tru_m);
    h3Jet_deltaR_MJ1_e3c_tru_m = new TH3F("h3Jet_deltaR_MJ1_e3c_tru_m", "h3Jet_deltaR_MJ1_e3c_tru_m",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3Jet_deltaR_MJ1_e3c_tru_m);
    h3Jet_deltaR_MJ2_e3c_tru_m = new TH3F("h3Jet_deltaR_MJ2_e3c_tru_m", "h3Jet_deltaR_MJ2_e3c_tru_m",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3Jet_deltaR_MJ2_e3c_tru_m);
    h3Jet_deltaR_MJ3_e3c_tru_m = new TH3F("h3Jet_deltaR_MJ3_e3c_tru_m", "h3Jet_deltaR_MJ3_e3c_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3Jet_deltaR_MJ3_e3c_tru_m);
    h3_JJMB_tru_m = new TH3F("h3_JJMB_tru_m", "h3_JJMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_JJMB_tru_m);
    h3_JMBMB_tru_m = new TH3F("h3_JMBMB_tru_m", "h3_JMBMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_JMBMB_tru_m);
    h3_MB1MB1MB2_tru_m = new TH3F("h3_MB1MB1MB2_tru_m", "h3_MB1MB1MB2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_MB1MB1MB2_tru_m);
    h3_MB1MB2MB2_tru_m = new TH3F("h3_MB1MB2MB2_tru_m", "h3_MB1MB2MB2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_MB1MB2MB2_tru_m);
    h3_JMB1MB2_tru_m = new TH3F("h3_JMB1MB2_tru_m", "h3_JMB1MB2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_JMB1MB2_tru_m);
    h3_MB1MB2MB3_tru_m = new TH3F("h3_MB1MB2MB3_tru_m", "h3_MB1MB2MB3_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_MB1MB2MB3_tru_m);
    h3_BMBMB_tru_m = new TH3F("h3_BMBMB_tru_m", "h3_BMBMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_BMBMB_tru_m);
    h3_SMBMB_tru_m = new TH3F("h3_SMBMB_tru_m", "h3_SMBMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_SMBMB_tru_m);
    h3_BMB1MB2_tru_m = new TH3F("h3_BMB1MB2_tru_m", "h3_BMB1MB2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_BMB1MB2_tru_m);
    h3_BBMB_tru_m = new TH3F("h3_BBMB_tru_m", "h3_BBMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_BBMB_tru_m);
    h3_SBMB_tru_m = new TH3F("h3_SBMB_tru_m", "h3_SBMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_SBMB_tru_m);
    h3_SSMB_tru_m = new TH3F("h3_SSMB_tru_m", "h3_SSMB_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_SSMB_tru_m);
    h3_SMB1MB2_tru_m = new TH3F("h3_SMB1MB2_tru_m", "h3_SMB1MB2_tru_m", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    // fOutput->Add(h3_SMB1MB2_tru_m);

if(fCout)cout<<"Unmatched jets only"<<endl;
    h3Jet_deltaR_MJ_e3c_um = new TH3F("h3Jet_deltaR_MJ_e3c_um", "h3Jet_deltaR_MJ_e3c_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ_e3c_um);
    h3Jet_deltaR_MJ0_e3c_um = new TH3F("h3Jet_deltaR_MJ0_e3c_um", "h3Jet_deltaR_MJ0_e3c_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ0_e3c_um);
    h3Jet_deltaR_MJ1_e3c_um = new TH3F("h3Jet_deltaR_MJ1_e3c_um", "h3Jet_deltaR_MJ1_e3c_um",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ1_e3c_um);
    h3Jet_deltaR_MJ2_e3c_um = new TH3F("h3Jet_deltaR_MJ2_e3c_um", "h3Jet_deltaR_MJ2_e3c_um",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ2_e3c_um);
    h3Jet_deltaR_MJ3_e3c_um = new TH3F("h3Jet_deltaR_MJ3_e3c_um", "h3Jet_deltaR_MJ3_e3c_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3Jet_deltaR_MJ3_e3c_um);
    
    h3_MB1MB1MB1_um = new TH3F("h3_MB1MB1MB1_um", "h3_MB1MB1MB1_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB1_um);
    h3_JJMB_um = new TH3F("h3_JJMB_um", "h3_JJMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JJMB_um);
    h3_JMBMB_um = new TH3F("h3_JMBMB_um", "h3_JMBMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMBMB_um);
    h3_MB1MB1MB2_um = new TH3F("h3_MB1MB1MB2_um", "h3_MB1MB1MB2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB1MB2_um);
    h3_MB1MB2MB2_um = new TH3F("h3_MB1MB2MB2_um", "h3_MB1MB2MB2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB2_um);
    h3_JMB1MB2_um = new TH3F("h3_JMB1MB2_um", "h3_JMB1MB2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_JMB1MB2_um);
    h3_MB1MB2MB3_um = new TH3F("h3_MB1MB2MB3_um", "h3_MB1MB2MB3_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_MB1MB2MB3_um);
    h3_BMBMB_um = new TH3F("h3_BMBMB_um", "h3_BMBMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BMBMB_um);
    h3_SMBMB_um = new TH3F("h3_SMBMB_um", "h3_SMBMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SMBMB_um);
    h3_BMB1MB2_um = new TH3F("h3_BMB1MB2_um", "h3_BMB1MB2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BMB1MB2_um);
    h3_SMB1MB2_um = new TH3F("h3_SMB1MB2_um", "h3_SMB1MB2_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SMB1MB2_um);
    h3_BBMB_um = new TH3F("h3_BBMB_um", "h3_BBMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_BBMB_um);
    h3_SBMB_um = new TH3F("h3_SBMB_um", "h3_SBMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SBMB_um);
    h3_SSMB_um = new TH3F("h3_SSMB_um", "h3_SSMB_um", nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(h3_SSMB_um);

  }

//   // /////////////////////////////////////////////////
   
    track_pt_matched = new TH2F("track_pt_match_hist", "Match Track Pt", 100, 0, 200, 100, 0, 200);
    // fOutput->Add(track_pt_matched);

    // Set track container pointers
    if(fCout) cout<<"Set track container pointers"<<endl;
    fDetectorLevel  = GetTrackContainer(fDetectorLevelName);
    fGeneratorLevel = GetMCParticleContainer(fGeneratorLevelName);
   
    //Unfolding optimization histograms from data
    if(!fDoEmbedding){
    if(fCout) cout<<"Unfolding optimization histograms from data"<<endl;
    OptUn_eec = new TH3F("Opt_Un_eec","Optimize unfolding bins",nJetPtbins, xbins, ndRbins, dRbins,nWtbins,wtbins);
    fOutput->Add(OptUn_eec);
    OptUn_e3c = new TH3F("Opt_Un_e3c","Optimize unfolding bins e3c",nJetPtbins, xbins, ndRbins, dRbins, nWtbins_e3c, wtbins_e3c);
    fOutput->Add(OptUn_e3c);
    }

    fTreeMatchTracks = new TTree("MatchTracksTree", "MatchTracksTree");
    if(fUnfolding==1)
     {
        if(fCout) cout<<"Unfolding tree"<<endl;
         
         fTreeMatchTracks->Branch("fJet_pt_det", &fJet_pt_det, "fJet_pt_det/F");
         fTreeMatchTracks->Branch("fJet_pt_tru", &fJet_pt_tru, "fJet_pt_tru/F");
         fTreeMatchTracks->Branch("fTrack_pt_det", &fTrack_pt_det, "fTrack_pt_det/F");
         fTreeMatchTracks->Branch("fTrack_pt_tru", &fTrack_pt_tru, "fTrack_pt_tru/F");
         fTreeMatchTracks->Branch("fTrack_eta_det", &fTrack_eta_det, "fTrack_eta_det/F");
         fTreeMatchTracks->Branch("fTrack_eta_tru", &fTrack_eta_tru, "fTrack_eta_tru/F");
         fTreeMatchTracks->Branch("fTrack_phi_det", &fTrack_phi_det, "fTrack_phi_det/F");
         fTreeMatchTracks->Branch("fTrack_phi_tru", &fTrack_phi_tru, "fTrack_phi_tru/F");

         fTreeMatchTracks->Branch("fTrack_pt_miss", &fTrack_pt_miss, "fTrack_pt_miss/F");
         fTreeMatchTracks->Branch("fTrack_eta_miss", &fTrack_eta_miss, "fTrack_eta_miss/F");
         fTreeMatchTracks->Branch("fTrack_phi_miss", &fTrack_phi_miss, "fTrack_phi_miss/F");

         fTreeMatchTracks->Branch("fTrack_pt_fake", &fTrack_pt_fake, "fTrack_pt_fake/F");
         fTreeMatchTracks->Branch("fTrack_eta_fake", &fTrack_eta_fake, "fTrack_eta_fake/F");
         fTreeMatchTracks->Branch("fTrack_phi_fake", &fTrack_phi_fake, "fTrack_phi_fake/F");
     }
    
    if(fCout) cout<<"histos and/or trees initialized"<<endl;

PostData(1, fOutput);
PostData(2, fTreeMatchTracks);  // Only post if the tree exists

}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().
  return kTRUE;
}

//_________________________________________________________________
//From AliAnalysisTaskPhiCorrelations
Long64_t AliAnalysisTaskJetsEECpbpb::GetUniqueEventID(AliVEvent* inputEvent)
{
  // Get event ID from header
  AliVHeader* eventIDHeader = inputEvent->GetHeader();
  if (eventIDHeader)
    return eventIDHeader->GetEventIdAsLong();
  else
    return 0;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::CountEmptyEvents()
{
  //
  // count Empty Events
  //
  Bool_t emptyEvent= kTRUE;
  UInt_t fAOD_FilterBits = 1<<8 | 1<<9;

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
    std::cout << " Info::Empty event " << std::endl;
  }
  else {
    if (fCout) std::cout << " Info:: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  }
  return emptyEvent;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::CheckHighPtTrackInEvent(AliVEvent* inputEvent)
{

  Bool_t skipEvent= kFALSE;
  UInt_t fAOD_FilterBits = 1<<8 | 1<<9;

  for (Int_t itrack=0;itrack<fAOD->GetNumberOfTracks();itrack++) {   // Track loop
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
    if(!track || !track->TestFilterBit(fAOD_FilterBits)) continue;
    //if a data track is greater than fHighPtTrackCutEvent, skip event
    if (track->Pt()> fHighPtTrackCutEvent && track->GetLabel()==-1) { 
      skipEvent= kTRUE;
      break;
    }
  }

  if (skipEvent) {
    if (fCout) std::cout << " Info::Event ignored due to high pT track " << std::endl;
  }
  else {
    if (fCout) std::cout << " Info:: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  }
  return skipEvent;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::FillHistograms()
{
//Embedding mode   
if (fDoEmbedding && fJetShapeType == kDetEmbPartPythia){
  Long64_t counter = GetUniqueEventID(InputEvent());
  if(fAddEventCuts)
  {
    if(CheckHighPtTrackInEvent(InputEvent()))
    {
      return 0;
    }
  }
  if(fCout){cout<<"Embedding loop"<<endl;}
  FillEmbHistograms(counter);
 }
else{
    
  AliEmcalJet *jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  // container zero is always the base containe: the data container, the
  // embedded subtracted in the case of embedding or the detector level in case
  // of pythia


  if (fCentSelectOn)
    if ((fCent > fCentMax) || (fCent < fCentMin))
      return 0;

  Float_t rhoVal = 0, rhoMassVal = 0.;
  if (jetCont) {
    jetCont->ResetCurrentID();
    if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kDerivSub)) {
      // rho
      AliRhoParameter *rhoParam = dynamic_cast<AliRhoParameter *>(
								  InputEvent()->FindListObject("RhoSparseR020"));
      if (!rhoParam) {
        Printf("%s: Could not retrieve rho %s (some histograms will be filled "
               "with zero)!",
               GetName(), jetCont->GetRhoName().Data());
      } else
        rhoVal = rhoParam->GetVal();
      // rhom
      AliRhoParameter *rhomParam = dynamic_cast<AliRhoParameter *>(
								   InputEvent()->FindListObject("RhoMassSparseR020"));
      if (!rhomParam) {
        Printf("%s: Could not retrieve rho_m %s (some histograms will be "
               "filled with zero)!",
               GetName(), jetCont->GetRhoMassName().Data());
      } else
        rhoMassVal = rhomParam->GetVal();
    }
    
    while ((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1)
        continue;
      AliEmcalJet *jet2 = 0x0;
      AliEmcalJet *jet3 = 0x0;
      AliEmcalJet *jetUnmatched = 0x0; //embedded level unmatched jet
      if (jet1->Pt() < fJetPtMin) {continue;} //Cuts on jet_pt
      AliEmcalJet *jetUS = NULL;
      Int_t ifound = 0, jfound = 0;
      Int_t ilab = -1, jlab = -1;

      int sub1 = -1;
      int sub2 = -1;

    //need this for pp analysis. Mode to run over pythia and match jets 
      if (fJetShapeType == kPythiaDef) {

        AliJetContainer *jetContTrue = GetJetContainer(1);
        AliJetContainer *jetContUS = GetJetContainer(2);
        AliJetContainer *jetContPart = GetJetContainer(3);
	
        if (fJetShapeSub == kConstSub) {
	  
          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1)
                ilab = i;
            }
          }
          if (ilab == -1)
            continue;
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();

          if (!jet2) {
            Printf("jet2 does not exist, returning");
            continue;
          }

          for (Int_t j = 0; j < jetContPart->GetNJets(); j++) {
	    
            jet3 = jetContPart->GetJet(j);
            if (!jet3)
              continue;
            if (jet3->GetLabel() == jet2->GetLabel()) {
              jfound++;
              if (jfound == 1)
                jlab = j;
            }
          }
          if (jlab == -1)
            continue;
          jet3 = jetContPart->GetJet(jlab);
          if (!jet3) {
            Printf("jet3 does not exist, returning");
            continue;
          }
        }
        if (!(fJetShapeSub == kConstSub)){
          jet3 = jet1->ClosestJet();
        }
        if (!jet3) {
          Printf("jet3 does not exist, returning");
          continue;
        }

        Double_t ptMatch=0;
        Int_t kMatched = 0;
            if (fJetShapeType == kPythiaDef)
            {
                kMatched = 1;
                if (fJetShapeSub == kConstSub)
                    kMatched = 3;
                ptMatch = jet3->Pt();
                //                cout<<"the matched jet "<<jet3->Pt()<<" "<<kMatched<<endl;
                // ComputeEncMC(jet1, jetCont, jet3, kMatched);
            }
    }
 
      Float_t ptSubtracted = 0;
      if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub) {  
        ptSubtracted = jet1->Pt();
     }
      else if (fJetShapeSub == kDerivSub) {
        ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
      }
    //we will almost exclusively use kNoSub
      else if (fJetShapeSub == kNoSub) {
        if ((fJetShapeType == kData))
          {
          ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
          if(ptSubtracted< fPtThreshold) continue;
          if(fCout) cout<<"data loop"<<endl;
          ComputeENC(jet1, ptSubtracted, jetCont);//Computing ENC on raw data
          if(fCout) cout<<"computed ENC in data loop"<<endl;
          FillDelPtCones(jet1,InputEvent(),GetRhoVal(0));
          std::vector<fastjet::PseudoJet> jetConstituents; 
          jetConstituents.clear();
          fastjet::PseudoJet PseudoTracksJet; //Creating a pseudojet object called PseduoTracks
          unsigned int constituentIndexJet = 0;
          for (auto partJet: jet1->GetParticleConstituents())
          {
          PseudoTracksJet.reset(partJet.Px(), partJet.Py(), partJet.Pz(), partJet.E());
          const AliVParticle* part2Jet = partJet.GetParticle(); 
          PseudoTracksJet.set_user_index(GetConstituentID(constituentIndexJet, part2Jet, jet1)); 
          if (PseudoTracksJet.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
          jetConstituents.push_back(PseudoTracksJet);
          constituentIndexJet++;
          }
          if(ifeec){
          if(fCout) cout<<"Now finding data cones for EEC"<<endl;
          auto cones = FindConesDataEEC(jet1,InputEvent());
          std::vector<fastjet::PseudoJet> cone1 = cones.first;
          std::vector<fastjet::PseudoJet> cone2 = cones.second;
          FillBkgSubJetsDataEEC(cone1,cone1,ptSubtracted,true,"sameMB");
          FillBkgSubJetsDataEEC(cone1,cone2,ptSubtracted,false,"diffMB");
          FillBkgSubJetsDataEEC(jetConstituents,cone1,ptSubtracted,false,"jetMB");
          }

          if(ife3c){
          if(fCout) cout<<"Now finding data cones for E3C"<<endl;
          std::vector<fastjet::PseudoJet> minBiasParticles, minBiasParticles2, minBiasParticles3;
          std::tie(minBiasParticles, minBiasParticles2, minBiasParticles3) = FindConesDataE3C(jet1,InputEvent());
          FillBkgSubJetsDataE3C(minBiasParticles,minBiasParticles,minBiasParticles,ptSubtracted,"all","sameMB");
          FillBkgSubJetsDataE3C(jetConstituents,jetConstituents,minBiasParticles,ptSubtracted,"two","jetjetMB");
          FillBkgSubJetsDataE3C(jetConstituents,minBiasParticles,minBiasParticles,ptSubtracted,"two","jetMBMB");
          FillBkgSubJetsDataE3C(minBiasParticles,minBiasParticles,minBiasParticles2,ptSubtracted,"two","MB1MB1MB2");
          FillBkgSubJetsDataE3C(minBiasParticles,minBiasParticles2,minBiasParticles2,ptSubtracted,"two","MB1MB2MB2");
          FillBkgSubJetsDataE3C(jetConstituents,minBiasParticles2,minBiasParticles3,ptSubtracted,"none","jetMB1MB2");
          FillBkgSubJetsDataE3C(minBiasParticles,minBiasParticles2,minBiasParticles3,ptSubtracted,"none","MB1MB2MB3");
          }
          }
        else if ((fJetShapeType == kPythiaDef) || (fJetShapeType == kMCTrue) || (fJetShapeType == kGenOnTheFly)){
          ptSubtracted = jet1->Pt();
      }
      }
    }
  } 
}
if(fCout){cout<<"Finished running things"<<endl;}

return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::GetTrueJetPtFraction(AliEmcalJet* jet, Double_t& truePtFraction, Double_t& truePtFraction_mcparticles)
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
  Double_t jetRadius = fConeR;
  if(fMCParticleArray)
  {
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
  }

  if(pt_all)
  {
    truePtFraction = (pt_truth/pt_all);
    truePtFraction_mcparticles = (pt_truth_mcparticles/pt_all);
  }
}

//________________________________________________________________________
//jet1 is hybrid jet, jet2 is matched truth jet
//Creates tree for unfolding - includes flags to ignore storing fake and missed tracks 
void AliAnalysisTaskJetsEECpbpb::FillMatchedTrackTree(AliEmcalJet* fJet_Hyb, AliEmcalJet* fJet_Tru, Float_t& JetEmbPtSub)
{
  if(fUnfolding==1)
  {
   //Hybrid level
    std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
    fConstituents.clear();
    fastjet::PseudoJet PseudoTracks;
    unsigned int JetconstituentIndex = 0;
    for (auto part: fJet_Hyb->GetParticleConstituents())
    {
        PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E());
        const AliVParticle* part2 = part.GetParticle(); 
        if(GetConstituentID(JetconstituentIndex, part2, fJet_Hyb)!= -1){
        PseudoTracks.set_user_index(GetConstituentID(JetconstituentIndex, part2, fJet_Hyb));
        JetconstituentIndex++;
        }
        else{
          PseudoTracks.set_user_index(-1); 
        }
        if (PseudoTracks.pt() < fMinENCtrackPt) continue; 
        fConstituents.push_back(PseudoTracks);
    }

    //Truth level 
    std::vector<fastjet::PseudoJet> fConstituents_tru; 
    fConstituents_tru.clear();
    fastjet::PseudoJet PseudoTracks_tru; 
    unsigned int JetconstituentIndex_tru = 0;
    for (auto part_tru: fJet_Tru->GetParticleConstituents())
    {
        PseudoTracks_tru.reset(part_tru.Px(), part_tru.Py(), part_tru.Pz(), part_tru.E()); 
        const AliVParticle* part_tru2 = part_tru.GetParticle(); 
        PseudoTracks_tru.set_user_index(GetConstituentID(JetconstituentIndex_tru, part_tru2, fJet_Tru));
        if(GetConstituentID(JetconstituentIndex_tru, part_tru2, fJet_Tru)!= -1){
        PseudoTracks_tru.set_user_index(GetConstituentID(JetconstituentIndex_tru, part_tru2, fJet_Tru));
        JetconstituentIndex_tru++;
        }
        else{
          PseudoTracks_tru.set_user_index(-1); 
        }

        if (PseudoTracks_tru.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents_tru.push_back(PseudoTracks_tru);
    }

    // double jet_pt = fJet_Hyb->Pt();
    float jet_pt = JetEmbPtSub; //need subtracted pt
    float jet_pt_tru = fJet_Tru->Pt();

    std::vector<int> tru_index,det_index;
    std::vector<fastjet::PseudoJet> matchtracks_det, matchtracks_tru, misstracks ,faketracks;
    matchtracks_tru.clear();
    matchtracks_det.clear();
    faketracks.clear();
    misstracks.clear();
    tru_index.clear();
    det_index.clear();

    for(int j=0; j<int(fConstituents.size()); j++)
    {
        det_index.push_back(fConstituents[j].user_index());
    }

    for(int i=0; i<int(fConstituents_tru.size()); i++)
    {
        tru_index.push_back(fConstituents_tru[i].user_index());
    }

    if(fFakeJetTrack==1)
    {
        for(int j = 0;j<int(fConstituents.size()); j++)
        {
            int valueToCheck = fConstituents[j].user_index();
            auto it = std::find(tru_index.begin(), tru_index.end(), valueToCheck);
            if (it != tru_index.end()){matchtracks_det.push_back(fConstituents[j]);}
            else {faketracks.push_back(fConstituents[j]);
                if(fUnfolding==1)
                {
                    fJet_pt_det = jet_pt;
                    fJet_pt_tru = jet_pt_tru;
                    fTrack_pt_tru = 0;
                    fTrack_eta_tru = 0;
                    fTrack_phi_tru = 0;
                    fTrack_pt_det = 0;
                    fTrack_eta_det = 0;
                    fTrack_phi_det = 0;
                    fTrack_eta_miss = 0;
                    fTrack_phi_miss = 0;
                    fTrack_pt_miss = 0;
                    fTrack_eta_fake = fConstituents[j].eta();
                    fTrack_phi_fake = fConstituents[j].phi();
                    fTrack_pt_fake = fConstituents[j].pt();
                    fTreeMatchTracks->Fill();
                }
            }
        }
    }

    if(fMissJetTrack==1)
    {
        for(int i=0; i<int(fConstituents_tru.size()); i++)
        {
            int valueToCheck = fConstituents_tru[i].user_index();
            auto it = std::find(det_index.begin(), det_index.end(), valueToCheck);
            if (it != det_index.end()){matchtracks_tru.push_back(fConstituents_tru[i]);}
            else {misstracks.push_back(fConstituents_tru[i]);
            if(fUnfolding==1)
                {
                    fJet_pt_tru = jet_pt_tru;
                    fJet_pt_det = jet_pt;
                    fTrack_pt_det = 0;
                    fTrack_eta_det = 0;
                    fTrack_phi_det = 0;
                    fTrack_eta_fake = 0;
                    fTrack_phi_fake = 0;
                    fTrack_pt_fake = 0;
                    fTrack_pt_tru = 0;
                    fTrack_eta_tru = 0;
                    fTrack_phi_tru = 0;
                    fTrack_eta_miss = fConstituents_tru[i].eta();
                    fTrack_phi_miss = fConstituents_tru[i].phi();
                    fTrack_pt_miss = fConstituents_tru[i].pt();
                    fTreeMatchTracks->Fill();
                }
            
            }
        }
    }

     //sort the tracks according to their indices
    struct CustomComparator {
        bool operator()(const fastjet::PseudoJet& track1, const fastjet::PseudoJet& track2) const {
            return track1.user_index() < track2.user_index();
        }
    };
  
    std::sort(matchtracks_tru.begin(), matchtracks_tru.end(), CustomComparator());
    std::sort(matchtracks_det.begin(), matchtracks_det.end(), CustomComparator());

    if(fMatchJetTrack==1)
    { 
      for(int j = 0; j < int(matchtracks_tru.size()); j++) // match
        {
         
                 fJet_pt_det = jet_pt;
                 fJet_pt_tru = jet_pt_tru;
                 fTrack_eta_miss = 0;
                 fTrack_phi_miss = 0;
                 fTrack_pt_miss = 0;
                 fTrack_eta_fake = 0;
                 fTrack_phi_fake = 0;
                 fTrack_pt_fake = 0;
                 fTrack_eta_tru = matchtracks_tru[j].eta();
                 fTrack_phi_tru = matchtracks_tru[j].phi();
                 fTrack_eta_det = matchtracks_det[j].eta();
                 fTrack_phi_det = matchtracks_det[j].phi();
                 fTrack_pt_det = matchtracks_det[j].pt();
                 fTrack_pt_tru = matchtracks_tru[j].pt();
                 fTreeMatchTracks->Fill();
             
        }
    }
  }
  else{
    if(fCout)cout<<"Not writing to trees because fUnfolding is set to FALSE"<<endl;
  }

  return;
}
//________________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::GetMatchedJetObservables(AliEmcalJet* jet, Double_t& detJetPt, Double_t& partJetPt, Double_t& detJetPhi, Double_t& detJetEta, Double_t& partJetPhi, Double_t& partJetEta, Double_t& detJetDistance, Double_t& partJetDistance, Float_t& JetEmbPtSub)
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

//Output trees 
if(fDoDetLevelMatching && fDoPartLevelMatching){
  if(fUnfolding == 1){FillMatchedTrackTree(jet,jet3,JetEmbPtSub);}
}
  return;
}
//
//________________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::FillEmbHistograms(Long64_t EvCounter)
{
 //
  Float_t pT_sub_min = fjetMinPtSub;
  if (fCout) cout << "Minimum jet pt after subtraction is " << fjetMinPtSub << endl;
  if (fCout) cout << "Maximum jet pt after subtraction is " << fjetMaxPtSub << endl;

  // Int_t trackOrigin = -10;

  // Float_t leadJetPhi_Emb = 999;
  // Float_t leadJetEta_Emb = 999;
  // Float_t leadJetPtSub_Emb = -999;

  // Float_t leadJetPhi_Rec = 999;
  // Float_t leadJetEta_Rec = 999;
  // Float_t leadJetPtSub_Rec = -999;

  // Float_t leadJetPhi_Gen = 999;
  // Float_t leadJetEta_Gen = 999;
  // Float_t leadJetPtSub_Gen = -999;

  //Gen Signal JETS
  // Get the jet container
  fGenJetContainer = this->GetJetContainer("GenJets"); //2 = Pythia particle truth jets
  if (this->GetJetContainer("GenJets")) fDoPartLevelMatching = kTRUE;
  TString fGenRhoName = fGenJetContainer->GetRhoName();
  if (fCout) cout << "Gen Rho Name is " << fGenRhoName << endl;

  if (fGenJetContainer->GetRhoParameter()) fjetGenRhoVal = fGenJetContainer->GetRhoVal();
  if (fCout) cout << "In the FillEmbJets Gen Rho value is " << fjetGenRhoVal << endl;

  //RECO Signal JETS
  // Get the jet container
  fRecoJetContainer = this->GetJetContainer("RecoJets"); //1 = Pythia det level jets
  if (this->GetJetContainer("RecoJets")) fDoDetLevelMatching = kTRUE; 
  TString fRecoRhoName = fRecoJetContainer->GetRhoName();
  if (fCout) cout << "Reco Rho Name is " << fRecoRhoName << endl;

  if (fRecoJetContainer->GetRhoParameter()) fjetRecoRhoVal = fRecoJetContainer->GetRhoVal();
  if (fCout) cout << "In the FillEmbJets Reco Rho value is " << fjetRecoRhoVal << endl;


  //EMB Signal JETS
  // Get the jet container
  fEmbJetContainer = this->GetJetContainer("EmbJets"); //0 = embedded = Data + Pythia det level jets
  TString fEmbRhoName = fEmbJetContainer->GetRhoName();
  if (fCout) cout << "Emb Rho Name is " << fEmbRhoName << endl;
  AliParticleContainer* embParticles = fEmbJetContainer->GetParticleContainer();
  if (embParticles) cout << "embParticles name is " << embParticles->GetName() << endl;


  if (fEmbJetContainer->GetRhoParameter()) fjetEmbRhoVal = fEmbJetContainer->GetRhoVal();
  if (fCout) cout << "In the FillEmbJets Emb Rho value is " << fjetEmbRhoVal << endl;

  if(fDoDetLevelMatching) DoJetMatching(); //this matches by R

  fEmbJetContainer->ResetCurrentID();

  for(auto jet_emb : fEmbJetContainer->accepted())
  {

    fhasAcceptedEMCjet = 1;
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

    //HERE is where the pT matching condition comes in
    GetMatchedJetObservables(jet_emb, matchedJetPt_Det, matchedJetPt_Part, matchedJetPhi_Det, matchedJetEta_Det, matchedJetPhi_Part, matchedJetEta_Part, matchedJetDistance_Det, matchedJetDistance_Part, jet_embptsub);
    
    std::vector<fastjet::PseudoJet> jetConstituents; 
    std::string matchedType = "";
    bool ifmj = false;
    if (matchedJetPt_Det==-0.1 || matchedJetPt_Part==-0.1) {
      if(jet_embptsub<fPtThreshold) continue;
      if(jet_embArea<fjetMinArea) continue;
      if(fCout){cout << "THIS JET WASNT MATCHED" << endl;} 
      hJet_num_um->Fill(0.5);
      h_jetpt_um->Fill(jet_embptsub);
      jetConstituents.clear();
      fastjet::PseudoJet PseudoJetTracks; 
      unsigned int JetconstituentIndex = 0;
      for (auto part: jet_emb->GetParticleConstituents())
      {
        PseudoJetTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); 
        const AliVParticle* part2 = part.GetParticle();
        //if its not a background particle, set ID 
        if(part2->GetLabel() == -1){
        PseudoJetTracks.set_user_index(-1);
        if(fCout)cout<<"THIS IS A bkg PARTICLE -- keeping index fixed to -1"<<endl;
        }
        else{
        PseudoJetTracks.set_user_index(GetConstituentID(JetconstituentIndex, part2, jet_emb));
        if(fCout)cout<<"THIS IS A signal PARTICLE -- setting index"<<endl;
        JetconstituentIndex++;
        }
        if (PseudoJetTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        jetConstituents.push_back(PseudoJetTracks);
       }
      if(jetConstituents.size() == 0) continue;
      matchedType = "unmatchedJet";
      if(fCout){cout<<"About to fill embedded jets"<<endl;}
        if(EvCounter % 2 == 0){
          if(fCout){cout<<"counter is even "<<EvCounter<<endl;}
           if(!ifcFactorHist){
            if(ifeec){
              if(fCout){cout<<"Unmatched Jet EEC cones"<<endl;}
            auto cones = FindThermalConeEEC(jet_emb,-2,-3);
            std::vector<fastjet::PseudoJet> minBiasParticles = cones.first;
            std::vector<fastjet::PseudoJet> minBiasParticles2 = cones.second;
            if(fCout){cout<<"Filling regular histograms "<<endl;}
            if(minBiasParticles.size()==0 || minBiasParticles2.size()==0){
              if(fCout)cout<<" cone 1 size "<< minBiasParticles.size()<<" and cone 2 size " << minBiasParticles2.size()<<endl;
              continue;
             } 
            FillEmbJetsEEC(jetConstituents, jetConstituents, jet_embptsub, jet_embptsub,true, "sameJet", -1, -1,false,matchedType.c_str());
            FillEmbJetsEEC(minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,true,"sameMB" , -2, -2,false,matchedType.c_str());
            FillEmbJetsEEC(minBiasParticles, minBiasParticles2, jet_embptsub, jet_embptsub,false,"diffMB", -2, -3,false,matchedType.c_str());
            FillEmbJetsEEC(jetConstituents, minBiasParticles, jet_embptsub, jet_embptsub,false,"jetMB", -1, -2,false,matchedType.c_str());
           }
           if(ife3c){
            if(fCout){cout<<"Unmatched Jet E3C cones"<<endl;}
            std::vector<fastjet::PseudoJet> minBiasParticles, minBiasParticles2, minBiasParticles3;
            std::tie(minBiasParticles, minBiasParticles2, minBiasParticles3) = FindThermalConeE3C(jet_emb,-2,-3,-4);
             if(minBiasParticles.size()==0 || minBiasParticles2.size()==0 || minBiasParticles3.size()==0 ){
              if(fCout)cout<<"one cone, 2 cones or all cones are empty, skipping "<<endl;
              continue;
             } 
            if(fCout){cout<<"######### Filling same jets E3C ######### "<<endl;}
            FillEmbJetsE3C(jetConstituents, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"all", "sameJet",ifmj);//jjj
            if(fCout){cout<<"######### Filling same mb : mb1mb1mb1 ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"all", "sameMB",ifmj);//mb^3
            if(fCout){cout<<"######### Filling two same : mb*jj  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb*j^2
            if(fCout){cout<<"######### Filling two same : mbmb*j  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents,minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//j*mb*mb
            if(fCout){cout<<"######### Filling two same : mb1mb2mb2  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2,minBiasParticles2, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb1mb2mb2
            if(fCout){cout<<"######### Filling two same : mb2mb2mb1  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles2, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb2mb2*mb1
            if(fCout){cout<<"######### Filling all diff : jmb1mb2  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents, minBiasParticles,minBiasParticles2, jet_embptsub, jet_embptsub,"alldiff", "",ifmj);//jet*mb1*mb2
            if(fCout){cout<<"######### Filling all diff : mb1mb2mb3  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2, minBiasParticles3, jet_embptsub, jet_embptsub,"alldiff", "", ifmj);
           }
          }
        }
        else{
        if(fCout){cout<<"counter is odd, in cFactor loop "<<EvCounter<<endl;}
        if(ifcFactorHist){
          if(ifeec){
            auto cones = FindThermalConeEEC(jet_emb,-2,-3);
            std::vector<fastjet::PseudoJet> minBiasParticles = cones.first;
            std::vector<fastjet::PseudoJet> minBiasParticles2 = cones.second;
          if(fCout){cout<<"counter is odd, cFactorHist is set to TRUE"<<endl;}
          if(minBiasParticles.size()==0 || minBiasParticles2.size()==0){
              if(fCout)cout<<" cone 1 size "<< minBiasParticles.size()<<" and cone 2 size " << minBiasParticles2.size()<<endl;
              continue;
             } 
          FillEmbJetsEEC(jetConstituents, jetConstituents, jet_embptsub, jet_embptsub,true, "sameJet", -1, -1,true,matchedType.c_str());
          FillEmbJetsEEC(minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,true,"sameMB" , -2, -2,true,matchedType.c_str());
          FillEmbJetsEEC(minBiasParticles, minBiasParticles2, jet_embptsub, jet_embptsub,false,"diffMB", -2, -3,true,matchedType.c_str());
          FillEmbJetsEEC(jetConstituents, minBiasParticles, jet_embptsub, jet_embptsub,false,"jetMB", -1, -2,true,matchedType.c_str());
          }
        if(ife3c){
             if(fCout){cout<<"Unmatched Jet E3C cones"<<endl;}
             std::vector<fastjet::PseudoJet> minBiasParticles, minBiasParticles2, minBiasParticles3;
             std::tie(minBiasParticles, minBiasParticles2, minBiasParticles3) = FindThermalConeE3C(jet_emb,-2,-3,-4);
              if(minBiasParticles.size()==0 || minBiasParticles2.size()==0 || minBiasParticles3.size()==0 ){
              if(fCout)cout<<"one cone, 2 cones or all cones are empty, skipping "<<endl;
              continue;
             } 
           if(fCout){cout<<"######### Filling same jets E3C ######### "<<endl;}
            FillEmbJetsE3C(jetConstituents, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"all", "sameJet",ifmj);//jjj
            if(fCout){cout<<"######### Filling same mb : mb1mb1mb1 ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"all", "sameMB",ifmj);//mb^3
            if(fCout){cout<<"######### Filling two same : mb*jj  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb*j^2
            if(fCout){cout<<"######### Filling two same : mbmb*j  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents,minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//j*mb*mb
            if(fCout){cout<<"######### Filling two same : mb1mb2mb2  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2,minBiasParticles2, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb1mb2mb2
            if(fCout){cout<<"######### Filling two same : mb2mb2mb1  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles2, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb2mb2*mb1
            if(fCout){cout<<"######### Filling all diff : jmb1mb2  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents, minBiasParticles,minBiasParticles2, jet_embptsub, jet_embptsub,"alldiff", "",ifmj);//jet*mb1*mb2
            if(fCout){cout<<"######### Filling all diff : mb1mb2mb3  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2, minBiasParticles3, jet_embptsub, jet_embptsub,"alldiff", "", ifmj);
           }
          }
        }
    }
    else{
      if(jet_embptsub<fPtThreshold) continue;
      if(jet_embArea<fjetMinArea) continue;
      float dpt = jet_embptsub - matchedJetPt_Det;
      h_dpt->Fill(dpt,jet_embptsub);
      h_jetpt_m->Fill(jet_embptsub);
      h_jetpt_m_tru->Fill(jet_embptsub,matchedJetPt_Part);
      if(fCout){cout << "THIS JET WAS MATCHED" << endl;} 
      jetConstituents.clear();
      fastjet::PseudoJet PseudoJetTracks; 
      unsigned int JetconstituentIndex = 0;
      for (auto part: jet_emb->GetParticleConstituents())
      {
        PseudoJetTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); 
        const AliVParticle* part2 = part.GetParticle(); 
        //if its not a background particle, set ID 
        if(part2->GetLabel() == -1){PseudoJetTracks.set_user_index(-1);}
        else{
        PseudoJetTracks.set_user_index(GetConstituentID(JetconstituentIndex, part2, jet_emb));
        JetconstituentIndex++;
        }
        if (PseudoJetTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        jetConstituents.push_back(PseudoJetTracks);
       }
       if(jetConstituents.size() == 0) continue;
       matchedType = "matchedJet";
       ifmj = true;
        if(EvCounter % 2 == 0){
          if(fCout){cout<<"counter is even "<<EvCounter<<endl;}
           if(!ifcFactorHist){
            if(ifeec){
              if(fCout){cout<<"MATCHED Jet EEC cones"<<endl;}
            auto cones = FindThermalConeEEC(jet_emb,-2,-3);
            std::vector<fastjet::PseudoJet> minBiasParticles = cones.first;
            std::vector<fastjet::PseudoJet> minBiasParticles2 = cones.second;
            if(fCout){cout<<"Filling regular histograms "<<endl;}
            if(minBiasParticles.size()==0 || minBiasParticles2.size()==0){
              if(fCout)cout<<" cone 1 size "<< minBiasParticles.size()<<" add cone 2 size " << minBiasParticles2.size()<<endl;
              continue;
             } 
            FillEmbJetsEEC(jetConstituents, jetConstituents, jet_embptsub, jet_embptsub,true, "sameJet", -1, -1,false,matchedType.c_str());
            FillEmbJetsEEC(minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,true,"sameMB" , -2, -2,false,matchedType.c_str());
            FillEmbJetsEEC(minBiasParticles, minBiasParticles2, jet_embptsub, jet_embptsub,false,"diffMB", -2, -3,false,matchedType.c_str());
            FillEmbJetsEEC(jetConstituents, minBiasParticles, jet_embptsub, jet_embptsub,false,"jetMB", -1, -2,false,matchedType.c_str());
           }
           if(ife3c){
            if(fCout){cout<<"MATCHED Jet E3C cones"<<endl;}
            std::vector<fastjet::PseudoJet> minBiasParticles, minBiasParticles2, minBiasParticles3;
            std::tie(minBiasParticles, minBiasParticles2, minBiasParticles3) = FindThermalConeE3C(jet_emb,-2,-3,-4);
            if(minBiasParticles.size()==0 || minBiasParticles2.size()==0 || minBiasParticles3.size()==0 ){
              if(fCout)cout<<"one cone, 2 cones or all cones are empty, skipping "<<endl;
              continue;
             } 
            if(fCout){cout<<"######### Filling same jets E3C ######### "<<endl;}
            FillEmbJetsE3C(jetConstituents, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"all", "sameJet",ifmj);//jjj
            if(fCout){cout<<"######### Filling same mb : mb1mb1mb1 ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"all", "sameMB",ifmj);//mb^3
            if(fCout){cout<<"######### Filling two same : mb*jj  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb*j^2
            if(fCout){cout<<"######### Filling two same : mbmb*j  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents,minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//j*mb*mb
            if(fCout){cout<<"######### Filling two same : mb1mb2mb2  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2,minBiasParticles2, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb1mb2mb2
            if(fCout){cout<<"######### Filling two same : mb2mb2mb1  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles2, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb2mb2*mb1
            if(fCout){cout<<"######### Filling all diff : jmb1mb2  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents, minBiasParticles,minBiasParticles2, jet_embptsub, jet_embptsub,"alldiff", "",ifmj);//jet*mb1*mb2
            if(fCout){cout<<"######### Filling all diff : mb1mb2mb3  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2, minBiasParticles3, jet_embptsub, jet_embptsub,"alldiff", "", ifmj);
           }
          }
        }
        else{
        if(fCout){cout<<"counter is odd, in cFactor loop "<<EvCounter<<endl;}
        if(ifcFactorHist){
          if(ifeec){
            auto cones = FindThermalConeEEC(jet_emb,-2,-3);
            std::vector<fastjet::PseudoJet> minBiasParticles = cones.first;
            std::vector<fastjet::PseudoJet> minBiasParticles2 = cones.second;
          if(minBiasParticles.size()==0 || minBiasParticles2.size()==0){
              if(fCout)cout<<" cone 1 size "<< minBiasParticles.size()<<" and cone 2 size " << minBiasParticles2.size()<<endl;
              continue;
             } 
          if(fCout){cout<<"counter is odd, cFactorHist is set to TRUE"<<endl;}
          FillEmbJetsEEC(jetConstituents, jetConstituents, jet_embptsub, jet_embptsub,true, "sameJet", -1, -1,true,matchedType.c_str());
          FillEmbJetsEEC(minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,true,"sameMB" , -2, -2,true,matchedType.c_str());
          FillEmbJetsEEC(minBiasParticles, minBiasParticles2, jet_embptsub, jet_embptsub,false,"diffMB", -2, -3,true,matchedType.c_str());
          FillEmbJetsEEC(jetConstituents, minBiasParticles, jet_embptsub, jet_embptsub,false,"jetMB", -1, -2,true,matchedType.c_str());
          }
        if(ife3c){
             if(fCout){cout<<"MATCHED Jet E3C cones"<<endl;}
             std::vector<fastjet::PseudoJet> minBiasParticles, minBiasParticles2, minBiasParticles3;
             std::tie(minBiasParticles, minBiasParticles2, minBiasParticles3) = FindThermalConeE3C(jet_emb,-2,-3,-4);
             if(minBiasParticles.size()==0 || minBiasParticles2.size()==0 || minBiasParticles3.size()==0 ){
              if(fCout)cout<<"one cone, 2 cones or all cones are empty, skipping "<<endl;
              continue;
             } 
           if(fCout){cout<<"######### Filling same jets E3C ######### "<<endl;}
            FillEmbJetsE3C(jetConstituents, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"all", "sameJet",ifmj);//jjj
            if(fCout){cout<<"######### Filling same mb : mb1mb1mb1 ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"all", "sameMB",ifmj);//mb^3
            if(fCout){cout<<"######### Filling two same : mb*jj  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, jetConstituents,jetConstituents, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb*j^2
            if(fCout){cout<<"######### Filling two same : mbmb*j  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents,minBiasParticles, minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//j*mb*mb
            if(fCout){cout<<"######### Filling two same : mb1mb2mb2  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2,minBiasParticles2, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb1mb2mb2
            if(fCout){cout<<"######### Filling two same : mb2mb2mb1  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles2, minBiasParticles,minBiasParticles, jet_embptsub, jet_embptsub,"two", "",ifmj);//mb2mb2*mb1
            if(fCout){cout<<"######### Filling all diff : jmb1mb2  ########### "<<endl;}
            FillEmbJetsE3C(jetConstituents, minBiasParticles,minBiasParticles2, jet_embptsub, jet_embptsub,"alldiff", "",ifmj);//jet*mb1*mb2
            if(fCout){cout<<"######### Filling all diff : mb1mb2mb3  ########### "<<endl;}
            FillEmbJetsE3C(minBiasParticles, minBiasParticles2, minBiasParticles3, jet_embptsub, jet_embptsub,"alldiff", "", ifmj);
           }
          }
        }
    }
  }
}
// //_______________________________________________________________________________________________
double AliAnalysisTaskJetsEECpbpb::delR(const fastjet::PseudoJet& ps1,const fastjet::PseudoJet& ps2)
    {
        double dphi = abs(ps1.phi()-ps2.phi());
        if(dphi>TMath::Pi()){dphi = (2.*TMath::Pi() - dphi);}
        double deta = ps1.eta()-ps2.eta();
        double dR = sqrt(dphi*dphi + deta*deta);
        return dR;
    }

// //_______________________________________________________________________________________________
double AliAnalysisTaskJetsEECpbpb::Calculate_pX( double pT, double eta, double phi)
{
	return(pT*TMath::Cos(phi));
}

double AliAnalysisTaskJetsEECpbpb::Calculate_pY( double pT, double eta, double phi)
{
	return(pT*TMath::Sin(phi));
}

double AliAnalysisTaskJetsEECpbpb::Calculate_pZ( double pT, double eta, double phi)
{
	return( pT*TMath::SinH(eta) );
}

double AliAnalysisTaskJetsEECpbpb::Calculate_E( double pT, double eta, double phi)
{
	double pZ = Calculate_pZ(pT,eta, phi);

	return( TMath::Sqrt(pT*pT + pZ*pZ) );
}

// //________________________________________________________________________
int AliAnalysisTaskJetsEECpbpb::GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet)
{
  // NOTE: Usually, we would use the global offset defined for the general subtracter extraction task. But we don't want to
  //       depend on that task, so we just define it here locally.
  //Get the id of the particle. If the label is not equal to -1, get the label and assign it to id. If it is -1, get the value of of track at id + 20000
  int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + 2000000);
  return id;
}


//_______________________________________________________________________
Double_t AliAnalysisTaskJetsEECpbpb::GetDownscaleWeight(string trigString)
{
  Double_t weight = 1.;
  TString triggerclass;
  if(trigString == "INT7") triggerclass = "CINT7-B-NOPF-CENT";
  else if(trigString == "EJ1") triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if(trigString == "EJ2") triggerclass = "CEMC7EJ2-B-NOPF-CENT";
  if(triggerclass.Length()) weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}

// //_______________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::FillDelPtCones(AliEmcalJet *fJet, AliVEvent* inputEvent, float rho){
 if (fCout) std::cout << "Creating del_pt cones for data jets" << std::endl;
    if (!fJet || !inputEvent) return; // Safety check

    double jetEta = fJet->Eta();
    double jet_phi = fJet->Phi();
    std::vector<fastjet::PseudoJet> tracksInCone;  
    Double_t Axis1 = 999, Axis2 = 999;
    float dPhi1, dPhi2, dEta, dPhi_plus, dPhi_minus;
    float trackptSum = 0.0; float trackptSum_bias1 = 0.0; float trackptSum_bias2 = 0.0;float ENCRandomConePt=0.0;
    float trackptSum_bias3 = 0.0; float trackptSum_bias4 = 0.0;
    bool fill_bias1 = false;
    bool fill_bias2 = false;
    bool fill_bias3 = false;
    bool fill_bias4 = false;
    if (fCout) std::cout << " Info:: ===== In the FillIncTracksReal ===== " << std::endl;

    UInt_t fAOD_FilterBits = (1 << 8) | (1 << 9);

    // Apply Random Cone or Perpendicular Cone logic
    if (fDoRandCone) {
        TRandom3 fRandom;
        double OffsetRandom = (1 + fRandom.Rndm()) * TMath::Pi() / 3; // Random Phi in range [π/3, 2π/3]
        Axis1 = jet_phi + OffsetRandom;
        Axis2 = jet_phi - OffsetRandom;
        if (Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi();
        if (Axis2 < 0) Axis2 += TMath::TwoPi();
    } 
    else if (fDoPerpCone) {
        Axis1 = ((jet_phi + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? jet_phi - ((3. / 2.) * TMath::Pi()) : jet_phi + (TMath::Pi() / 2.);
    }

    // Get the event
    AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(inputEvent);
    if (!fAOD) return;
    
    if (fCout) std::cout << "Number of tracks in event: " << fAOD->GetNumberOfTracks() << std::endl;

    // Track loop
    for (Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); ++itrack) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
        if (!track) continue;

        // Apply filter bit cuts
        if (!track->TestFilterBit(fAOD_FilterBits)) continue;

        // Apply track pT cuts
        if (track->Pt() > fMaxPtTrack) continue;
        if (fCout) std::cout << "Track passed pT cuts" << std::endl;
        
        Float_t track_phi = track->Phi();
        Float_t track_eta = track->Eta();

        // Compute dPhi and dEta
        dPhi1 = TMath::Abs(track_phi - Axis1);
        dPhi1 = (dPhi1 > TMath::Pi()) ? (2 * TMath::Pi() - dPhi1) : dPhi1;
        dEta = jetEta - track->Eta();

        if ((dPhi1 * dPhi1 + dEta * dEta) < (fConeR*fConeR)){

          if (track->Pt() > fMinENCtrackPt) {fill_bias1 = true;}//cone has one track greater than fMinENCtrackPt
          if (track->Pt() > fMinPtConst){fill_bias2 = true;}//cone has one track greater than fMinPtConst
          if (track->Pt() > fMinPtConst-2){fill_bias3 = true;}//cone has one track greater than fMinPtConst - 2
          if (track->Pt() > fMinPtConst+2){fill_bias4 = true;}//cone has one track greater than fMinPtConst + 2

         // Always add track to total sum
            trackptSum += track->Pt();

            if (track->Pt() > fMinENCtrackPt) ENCRandomConePt += track->Pt();
        }
    }

    // If a flag is set, sum all tracks for that bias
    if (fill_bias1) trackptSum_bias1 = trackptSum;
    if (fill_bias2) trackptSum_bias2 = trackptSum;
    if (fill_bias3) trackptSum_bias3 = trackptSum;
    if (fill_bias4) trackptSum_bias4 = trackptSum;

    float rhoArea = rho*(TMath::Pi())*fConeR*fConeR;

    delta_pt_cone->Fill(trackptSum - rhoArea);
    if(fill_bias1){
    delta_pt_coneBias1->Fill(trackptSum_bias1 - rhoArea);
    }
    if(fill_bias2){
    delta_pt_coneBias2->Fill(trackptSum_bias2 - rhoArea);
    }
    if(fill_bias3){
    delta_pt_coneBias3->Fill(trackptSum_bias3 - rhoArea);
    }
    if(fill_bias4){
    delta_pt_coneBias4->Fill(trackptSum_bias4 - rhoArea);
    }

    delta_pt_ENCcone->Fill(ENCRandomConePt - rhoArea);

}
// //_______________________________________________________________________
// //Perform jet matching on embedded jets
void AliAnalysisTaskJetsEECpbpb::DoJetMatching(){
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
  PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fMatchR);
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
  PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fMatchR);
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
  PerformGeometricalJetMatching(*jetsDetLevel, *jetsPartLevel, fMatchR);
}

//_______________________________________________________________________
//Use geometric matching for embedded jets
bool AliAnalysisTaskJetsEECpbpb::PerformGeometricalJetMatching(AliJetContainer& contBase,
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
//______________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::MatchTracks()
{
    AliVTrack* track;
    for (auto trackIterator : fDetectorLevel->accepted_momentum() )
    {
        track = trackIterator.second;
        Byte_t type = fDetectorLevel->GetTrackType(track);
        if (type <= 2) {
            Double_t sigma = 0;

            if (fIsEsd) {
                continue;
            } else { // AOD

                AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
                if(!aodtrack) AliFatal("Not a standard AOD");

                Int_t label = TMath::Abs(track->GetLabel());

                //            FillDetectorLevelTHnSparse(fCent, track->Eta(), track->Phi(), track->Pt(), sigma, type);

                if (fGeneratorLevel && label > 0) {
                    AliAODMCParticle *part =  fGeneratorLevel->GetAcceptMCParticleWithLabel(label); //label is used to match the generator to detector level
                    if (part) {
                        if (part->GetGeneratorIndex() == 0) {
                                track_pt_matched->Fill(part->Pt(),track->Pt());
                        }
                    }
                }
            }
        }
        else {
            AliError(Form("Track %d has type %d not recognized!", fDetectorLevel->GetCurrentID(), type));
        }
    }
}
//______________________________________________________________________
//To perform background subtraction in data 
//This will give BB' and JB contributions
void AliAnalysisTaskJetsEECpbpb::FillBkgSubJetsDataEEC(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, float jetpt, bool typeSame, std::string type)
{
  float w_eec = 0.0;
  float w_eec_3D = 0.0;

    if(typeSame){
        int mult = particles.size();
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = i+1; j < mult; j++)
            {
                if(particles.at(j).pt()<fMinENCtrackPt) continue;

                  float pt_i = particles.at(i).pt();
                  float pt_j =  particles.at(j).pt();
                
                  w_eec = (1./(jetpt*jetpt))*(2 * pt_i * pt_j);
                  w_eec_3D =  (1.0 /(jetpt*jetpt))*(pt_i * pt_j);

                  float dR = delR(particles.at(i), particles.at(j));

                  if(type == "sameMB"){
                    h_MB1_dat->Fill(jetpt, dR, w_eec);
                   
                    h3Jet_deltaR_MB1_dat->Fill(jetpt, dR, w_eec_3D);
                    h3Jet_deltaR_MB1_dat->Fill(jetpt, dR, w_eec_3D);
                  }

            }    
        }
    }
    else{
       int mult = particles.size();
       int mult2 = particles2.size();
      for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = 0; j < mult2; j++)
            {
                 if(particles2.at(j).pt()<fMinENCtrackPt) continue;
                  
                  float pt_i = particles.at(i).pt();
                  float pt_j =  particles2.at(j).pt();
                    
                  w_eec = (1./(jetpt*jetpt))*(2*pt_i * pt_j);
                  w_eec_3D =  (1.0 /(jetpt*jetpt))*(pt_i * pt_j);
                
                  float dR = delR(particles.at(i), particles2.at(j));

                if(type == "jetMB"){
                    h_JMB_dat->Fill(jetpt, dR, w_eec);

                    h3Jet_deltaR_JMB_dat->Fill(jetpt, dR, w_eec_3D);
                    h3Jet_deltaR_JMB_dat->Fill(jetpt, dR, w_eec_3D);
                }  
                if(type == "diffMB"){
                    h_MB1MB2_dat->Fill(jetpt, dR, w_eec);
                    
                    h3Jet_deltaR_MB1MB2_dat->Fill(jetpt, dR, w_eec_3D);
                    h3Jet_deltaR_MB1MB2_dat->Fill(jetpt, dR, w_eec_3D); 

                  }
                }
            }
        }   
  }
//________________________________________________________________________
//To perform background subtraction in data E3C
//This will give B'B'B' (all same)
//B'B''B'', B'B'B'' and JB'B',  JJB' (two diff)
//JB'B'', B'B''B''' (all diff)
void AliAnalysisTaskJetsEECpbpb::FillBkgSubJetsDataE3C(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2,std::vector<fastjet::PseudoJet> particles3, float jetpt, std::string typeSame, std::string type)
{
  int mult = particles.size();
  int mult2 = particles2.size();
  int mult3 = particles3.size();
    
  float w = 0;
  float w_iij = 0;
  float w_e3c_3D = 0;
  float w_iij_3D = 0;
    
  float w_ijs = 0;
  float w_ijs_3D = 0;
    
  float w_twosame = 0;
  float w_twosame_3D = 0;
    
  float w_ijk = 0;
  float w_ijk_3D = 0;

    if(typeSame == "all"){
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = i+1; j < mult; j++)
            {
                if(particles.at(j).pt()<fMinENCtrackPt) continue;
                
                if(fCout){cout<<"in this loop"<<endl;}
                
                w = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt());
                w_iij = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt());
                w_e3c_3D =  (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt());
                w_iij_3D =  (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt());

                float dR = delR(particles.at(i), particles.at(j));

                 if(type == "sameMB"){
                    h_MB1MB1MB1_dat->Fill(jetpt, dR, w);
                    h_MB1MB1MB1_dat->Fill(jetpt, dR, w_iij);
                    
                    h3_MB1MB1MB1_dat->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_dat->Fill(jetpt, dR, w_iij_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_dat->Fill(jetpt, dR, w_iij_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_dat->Fill(jetpt, dR, w_iij_3D);
                 }

                 for (int s = j+1; s < mult; s++)
                {
                    
                    if(particles.at(s).pt()<fMinENCtrackPt) continue;
                    w_ijs = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt());
                    w_ijs_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt());

                    float dR_ij = delR(particles.at(i), particles.at(j));
                    float dR_js = delR(particles.at(j), particles.at(s));
                    float dR_is = delR(particles.at(s), particles.at(i));
                    
                    float R_L = -1;
                    if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
                    else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
                    else{R_L = dR_is;}
                    
                    h_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs);
                            
                    h3_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs_3D);
                    h3_MB1MB1MB1_dat->Fill(jetpt, R_L, w_ijs_3D);

                }
            }
        }
    }  
    else if (typeSame == "two")
    {
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int  j = 0; j < mult2; j++)
            {
                
              if(particles2.at(j).pt()<fMinENCtrackPt) continue;

              w_twosame = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt());
              w_twosame_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt());
             
              float R_L = delR(particles.at(i), particles2.at(j));

              if(type == "jetjetMB"){
                h_JJMB_dat->Fill(jetpt, R_L,w_twosame);
                    
                h3_JJMB_dat->Fill(jetpt, R_L,w_twosame_3D);
                h3_JJMB_dat->Fill(jetpt, R_L,w_twosame_3D);
                h3_JJMB_dat->Fill(jetpt, R_L,w_twosame_3D);

              }
              else if(type == "jetMBMB"){
                h_JMBMB_dat->Fill(jetpt, R_L,w_twosame);
                    
                h3_JMBMB_dat->Fill(jetpt, R_L,w_twosame_3D);
                h3_JMBMB_dat->Fill(jetpt, R_L,w_twosame_3D);
                h3_JMBMB_dat->Fill(jetpt, R_L,w_twosame_3D);

              }
              else if(type == "MB1MB1MB2"){
                h_MB1MB1MB2_dat->Fill(jetpt, R_L, w_twosame);
                    
                h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_twosame_3D);
                h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_twosame_3D);
                h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_twosame_3D);

              }
              else {
                h_MB1MB2MB2_dat->Fill(jetpt, R_L, w_twosame);
                    
                h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_twosame_3D);
                h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_twosame_3D);
                h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_twosame_3D);
              }
                for (int k = j+1; k < mult2; k++)
                   {
                    if(particles2.at(k).pt()<fMinENCtrackPt) continue;
                    w_ijk = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt());
                    w_ijk_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt());

                    float dR_ij = delR(particles.at(i), particles2.at(j));
                    float dR_js = delR(particles2.at(j), particles2.at(k));
                    float dR_is = delR(particles2.at(k), particles.at(i));
                    
                    float R_L = -1;
                    if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
                    else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
                    else{R_L = dR_is;}

                    if(type == "jetjetMB"){
                      h_JJMB_dat->Fill(jetpt, R_L,w_ijk);
                          
                      h3_JJMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JJMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JJMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JJMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JJMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JJMB_dat->Fill(jetpt, R_L,w_ijk_3D);

                    }
                    else if(type == "jetMBMB"){
                      h_JMBMB_dat->Fill(jetpt, R_L,w_ijk);
                          
                      h3_JMBMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMBMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMBMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMBMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMBMB_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMBMB_dat->Fill(jetpt, R_L,w_ijk_3D);

                    }
                    else if(type == "MB1MB1MB2"){
                      h_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk);
                          
                      h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB1MB2_dat->Fill(jetpt, R_L, w_ijk_3D);

                    }
                    else {
                      h_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk);
                          
                      h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                      h3_MB1MB2MB2_dat->Fill(jetpt, R_L, w_ijk_3D);
                    }

                
                  }

            }    
        }
    }
    else{
       for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
           
            for (int  j = 0; j < mult2; j++)
            {
                if(particles2.at(j).pt()<fMinENCtrackPt) continue;
               
                for (int k = 0; k < mult3; k++)
                {
                    if(particles3.at(k).pt()<fMinENCtrackPt) continue;
                    
                    w_ijk = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt());
                    w_ijk_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt());
                    
                    float dR_ij = delR(particles.at(i), particles2.at(j));
                    float dR_js = delR(particles2.at(j), particles3.at(k));
                    float dR_is = delR(particles.at(i), particles3.at(k));
                    
                    float R_L = -1;
                    if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
                    else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
                    else{R_L = dR_is;}

                    if(type == "jetMB1MB2"){
                      h_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk);
                        
                      h3_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_JMB1MB2_dat->Fill(jetpt, R_L,w_ijk_3D);

                    }
                    if(type == "MB1MB2MB3"){
                      h_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk);
                        
                      h3_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk_3D);
                      h3_MB1MB2MB3_dat->Fill(jetpt, R_L,w_ijk_3D);

                    }
                }
            }
        }
    }
}
//________________________________________________________________________
// //Find cones in data to perform background subtraction with 
std::pair<std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>> 
AliAnalysisTaskJetsEECpbpb::FindConesDataEEC(AliEmcalJet *fJet, AliVEvent* inputEvent) {
    if (fCout) std::cout << "Creating cones for EEC data jets" << std::endl;
    if (!fJet || !inputEvent) return {}; // Safety check

    double jetEta = fJet->Eta();
    double jet_phi = fJet->Phi();

    Double_t Axis1 = 999, Axis2 = 999;
    float dPhi1, dPhi2, dEta, dPhi_plus, dPhi_minus;
    
    std::vector<fastjet::PseudoJet> tracksInConePlus;
    std::vector<fastjet::PseudoJet> tracksInConeMinus;

    if (fCout) std::cout << " Info:: ===== In the FillIncTracksReal ===== " << std::endl;

    UInt_t fAOD_FilterBits = (1 << 8) | (1 << 9);

    // Apply Random Cone or Perpendicular Cone logic
    if (fDoRandCone) {
        TRandom3 fRandom;
        double OffsetRandom = (1 + fRandom.Rndm()) * TMath::Pi() / 3; // Random Phi in range [π/3, 2π/3]
        Axis1 = jet_phi + OffsetRandom;
        Axis2 = jet_phi - OffsetRandom;
        if (Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi();
        if (Axis2 < 0) Axis2 += TMath::TwoPi();
    } 
    else if (fDoPerpCone) {
        Axis1 = ((jet_phi + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? jet_phi - ((3. / 2.) * TMath::Pi()) : jet_phi + (TMath::Pi() / 2.);
        Axis2 = ((jet_phi - (TMath::Pi() / 2.)) < 0) ? jet_phi + ((3. / 2.) * TMath::Pi()) : jet_phi - (TMath::Pi() / 2.);
    }

    // Get the event
    AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(inputEvent);
    if (!fAOD) return {};
    
    if (fCout) std::cout << "Number of tracks in event: " << fAOD->GetNumberOfTracks() << std::endl;

    // Track loop
    for (Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); ++itrack) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
        if (!track) continue;

        // Apply filter bit cuts
        if (!track->TestFilterBit(fAOD_FilterBits)) continue;

        // Apply track pT cuts
        if (track->Pt() < fMinENCtrackPt || track->Pt() > fMaxPtTrack) continue;
        if (fCout) std::cout << "Track passed pT cuts" << std::endl;
        
        Float_t track_phi = track->Phi();
        Float_t track_eta = track->Eta();

        // Compute dPhi and dEta
        dPhi1 = TMath::Abs(track_phi - Axis1);
        dPhi1 = (dPhi1 > TMath::Pi()) ? (2 * TMath::Pi() - dPhi1) : dPhi1;
        dPhi2 = TMath::Abs(track_phi - Axis2);
        dPhi2 = (dPhi2 > TMath::Pi()) ? (2 * TMath::Pi() - dPhi2) : dPhi2;

        dEta = jetEta - track->Eta();

        if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) < fConeR)) {
            Float_t rotated_phi = jet_phi + (track_phi - Axis1);
            if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
            if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
            fastjet::PseudoJet pseudoTrack(track->Pt() * TMath::Cos(rotated_phi), 
                                                track->Pt() * TMath::Sin(rotated_phi),
                                                track->Pt() * TMath::SinH(track_eta),
                                                track->E());
            tracksInConePlus.push_back(pseudoTrack);
          }

        if ((TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) < fConeR)) {
              Float_t rotated_phi = jet_phi + (track_phi - Axis2);
              if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
              if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();

              fastjet::PseudoJet pseudoTrack(track->Pt() * TMath::Cos(rotated_phi), 
                                               track->Pt() * TMath::Sin(rotated_phi),
                                               track->Pt() * TMath::SinH(track_eta),
                                               track->E());
              tracksInConeMinus.push_back(pseudoTrack);
        }

        if (fCout) std::cout << "Tracks in cone sizes: Plus = " 
                             << tracksInConePlus.size() 
                             << ", Minus = " 
                             << tracksInConeMinus.size() << std::endl;
    }

    return {tracksInConePlus, tracksInConeMinus}; // Return both vectors as a pair
}

//________________________________________________________________________
std::tuple<std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>>
AliAnalysisTaskJetsEECpbpb::FindConesDataE3C(AliEmcalJet *fJet, AliVEvent* inputEvent)
{
    if(fCout) cout<<"Creating cones for E3C data jets" <<endl;
    std::vector<fastjet::PseudoJet> coneParticles1, coneParticles2, coneParticles3; // Separate vectors for each cone

    // Jet kinematics
    Float_t jet_phi = fJet->Phi();
    Float_t jet_eta = fJet->Eta();

    // Anti-jet location (π - jet_phi)
    Float_t anti_jet_phi = jet_phi + TMath::Pi();
    if (anti_jet_phi > TMath::TwoPi()) anti_jet_phi -= TMath::TwoPi();

    Double_t Axis1_Perp, Axis2_Shifted, Axis3_Offset;
    Double_t dPhi1, dPhi2, dPhi3, dEta1, dEta2, dEta3;

    // **Cone 1: Perpendicular to jet axis (π/2 shift)**
    Axis1_Perp = jet_phi + (TMath::Pi() / 2.);
    if (Axis1_Perp > TMath::TwoPi()) Axis1_Perp -= TMath::TwoPi();

    // **Cone 2: Shift φ by fDeltaAxisShift and reflect η relative to the jet**
    Axis2_Shifted = jet_phi + fDeltaAxisShift;
    if (Axis2_Shifted > TMath::TwoPi()) Axis2_Shifted -= TMath::TwoPi();
    Float_t eta_reflected = -jet_eta; // Use reflected η for track selection

    // **Cone 3: fDeltaAxisShift away from Cone 2, same η as the jet**
    Axis3_Offset = Axis2_Shifted + fDeltaAxisShift;
    if (Axis3_Offset > TMath::TwoPi()) Axis3_Offset -= TMath::TwoPi();

    if (fCout) std::cout << " Info:: ===== In the FillIncTracksReal for E3C ===== " << std::endl;

    UInt_t fAOD_FilterBits = 1<<8 | 1<<9;

    // Get the event
    AliVEvent *event = inputEvent;
    if (!event) return {};
    fisGoodIncEvent = 1;
    fAOD = dynamic_cast<AliAODEvent*>(event);

    // --------------------------------------------------------------
    //  Main track loop
    // --------------------------------------------------------------

    for (Int_t itrack = 0; itrack < event->GetNumberOfTracks(); ++itrack) {   // Track loop
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
        if (!track) continue;

        // Apply filter cuts
        if (!track->TestFilterBit(fAOD_FilterBits)) continue;

        // Apply pT cuts
        if (track->Pt() < fMinENCtrackPt || track->Pt() > fMaxPtTrack) continue;

        // Store original track η and φ for selection
        Float_t track_phi = track->Phi();
        Float_t track_eta = track->Eta();

        // Compute distances to cones using original η-φ values
        dPhi1 = TMath::Abs(track_phi - Axis1_Perp);
        dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
        dEta1 = jet_eta - track_eta;
        double distanceCone1 = TMath::Sqrt(dPhi1 * dPhi1 + dEta1 * dEta1);

        dPhi2 = TMath::Abs(track_phi - Axis2_Shifted);
        dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
        dEta2 = eta_reflected - track_eta;
        double distanceCone2 = TMath::Sqrt(dPhi2 * dPhi2 + dEta2 * dEta2);

        dPhi3 = TMath::Abs(track_phi - Axis3_Offset);
        dPhi3 = (dPhi3 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi3 : dPhi3;
        dEta3 = jet_eta - track_eta;
        double distanceCone3 = TMath::Sqrt(dPhi3 * dPhi3 + dEta3 * dEta3);


        // Assign to correct cone
        if (distanceCone1 < fConeR) {
          Float_t rotated_phi = jet_phi + (track_phi - Axis1_Perp);  // Rotation in φ
          if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
          if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
          fastjet::PseudoJet PseudoTracks(track->Pt() * TMath::Cos(rotated_phi), 
                                          track->Pt() * TMath::Sin(rotated_phi),
                                          track->Pt() * TMath::SinH(track_eta), // Keeping η the same
                                          track->E());
          coneParticles1.push_back(PseudoTracks);
          }
        if (distanceCone2 < fConeR) {
          Float_t rotated_phi = jet_phi + (track_phi - Axis2_Shifted);  // Rotation in φ
          if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
          if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
          Float_t new_eta = -track_eta;
          fastjet::PseudoJet PseudoTracks(track->Pt() * TMath::Cos(rotated_phi), 
                                          track->Pt() * TMath::Sin(rotated_phi),
                                          track->Pt() * TMath::SinH(new_eta), 
                                          track->E());     
          coneParticles2.push_back(PseudoTracks);
          }
        if (distanceCone3 < fConeR) {
          Float_t rotated_phi = jet_phi + (track_phi - Axis3_Offset);  // Rotation in φ
          if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
          if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
          fastjet::PseudoJet PseudoTracks(track->Pt() * TMath::Cos(rotated_phi), 
                                          track->Pt() * TMath::Sin(rotated_phi),
                                          track->Pt() * TMath::SinH(track_eta), // Keeping η the same
                                          track->E());
          coneParticles3.push_back(PseudoTracks);
        }
    }

    return std::make_tuple(coneParticles1, coneParticles2, coneParticles3);
}



//________________________________________________________________________
//bkgindex = -1 for jet background particle, -2 for a thermal cone, -3 for second thermal cone
//pt is true (gen) level jet pT, jetpt is the pT of the embedded subtracted jet
//if cfactor && ifcFactorHist is true then fill the c-factor histograms
//this will be true 50% of times and false 50% times and will be set in the embedding
//if you want to use all the same events to find the c-factor just set it manually
void AliAnalysisTaskJetsEECpbpb::FillEmbJetsEEC(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, float jetpt, float pt, bool typeSame, std::string type, double bkgIndex1, double bkgIndex2, bool cfactor, std::string tag)
{
float w_eec = 0.0;
float w_eec_tru = 0.0;

float w_eec_3D = 0.0;
float w_eec_tru_3D = 0.0;

if(fCout) cout<<"In FillEmbJetsEEC"<<endl;
    if(typeSame)
    {
        int mult = particles.size();
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = i+1; j < mult; j++)
            {
                if(particles.at(j).pt()<fMinENCtrackPt) continue;

                  float pt_i = particles.at(i).pt();
                  float pt_j =  particles.at(j).pt();
                
                  w_eec = (1./(jetpt*jetpt))*(2 * pt_i * pt_j);
                  w_eec_tru = (1./(pt*pt))*(2 * pt_i * pt_j);

                  w_eec_3D =  (1.0 /(jetpt*jetpt))*(pt_i * pt_j);
                  w_eec_tru_3D =  (1.0 / (pt * pt))* (pt_i * pt_j);
                
                
                float dR = delR(particles.at(i), particles.at(j));
                int i1 = particles.at(i).user_index();
                int i2 = particles.at(j).user_index();
                
                if(type == "sameJet")
                {
                  if(fCout)std::cout<<"In type: sameJet"<<endl;
                        h_MJ->Fill(jetpt, pt, dR, w_eec);

                        h3Jet_deltaR_MJ->Fill(jetpt, dR, w_eec_3D);
                        h3Jet_deltaR_MJ->Fill(jetpt, dR, w_eec_3D);

                        if(tag=="matchedJet")
                        {
                          h_MJ_m->Fill(jetpt, pt, dR, w_eec);
                          h_MJ_tru_m->Fill(jetpt, pt, dR, w_eec_tru); //with true jet pT
                         
                          h3Jet_deltaR_MJ_m->Fill(jetpt, dR, w_eec_3D);
                          h3Jet_deltaR_MJ_m->Fill(jetpt, dR, w_eec_3D);

                          h3Jet_deltaR_MJ_tru_m->Fill(pt, dR, w_eec_tru_3D);
                          h3Jet_deltaR_MJ_tru_m->Fill(pt, dR, w_eec_tru_3D);

                          if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_MJ_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ_minpt->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MJ_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MJ_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                          }
                        }
                        else
                        {
                          h_MJ_um->Fill(jetpt, pt, dR, w_eec);
                        
                          h3Jet_deltaR_MJ_um->Fill(jetpt, dR, w_eec_3D);
                          h3Jet_deltaR_MJ_um->Fill(jetpt, dR, w_eec_3D);
                        }
                    
                    if ((i1 == bkgIndex1) && (i2 == bkgIndex1)) 
                    { 
                          h_MJ0->Fill(jetpt, pt, dR, w_eec);//both background

                          h3Jet_deltaR_MJ0->Fill(jetpt, dR, w_eec_3D);
                          h3Jet_deltaR_MJ0->Fill(jetpt, dR, w_eec_3D);

                          if(tag=="matchedJet"){
                            h_MJ0_m->Fill(jetpt, pt, dR, w_eec);//both background
                            h_MJ0_tru_m->Fill(jetpt, pt, dR, w_eec_tru);//both background
                            h3Jet_deltaR_MJ0_m->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ0_m->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MJ0_tru_m->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MJ0_tru_m->Fill(pt, dR, w_eec_tru_3D);

                              if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                                h3Jet_deltaR_MJ0_minpt->Fill(jetpt, dR, w_eec_3D);
                                h3Jet_deltaR_MJ0_minpt->Fill(jetpt, dR, w_eec_3D);

                                h3Jet_deltaR_MJ0_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                                h3Jet_deltaR_MJ0_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                              }
                          }
                          else{
                            h_MJ0_um->Fill(jetpt, pt, dR, w_eec);//both background
                            h3Jet_deltaR_MJ0_um->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ0_um->Fill(jetpt, dR, w_eec_3D);
                          }
                        }
                    else if ((i1 == bkgIndex1) || (i2 == bkgIndex1)){
                   
                          h_MJ1->Fill(jetpt, pt, dR, w_eec);// one is pythia

                          h3Jet_deltaR_MJ1->Fill(jetpt, dR, w_eec_3D);
                          h3Jet_deltaR_MJ1->Fill(jetpt, dR, w_eec_3D);

                          if(tag=="matchedJet")
                          {
                            h_MJ1_m->Fill(jetpt, pt, dR, w_eec);// one is pythia
                            h_MJ1_tru_m->Fill(jetpt, pt, dR, w_eec_tru);// one is pythia

                            h3Jet_deltaR_MJ1_m->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ1_m->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MJ1_tru_m->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MJ1_tru_m->Fill(pt, dR, w_eec_tru_3D);

                            if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_MJ1_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ1_minpt->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MJ1_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MJ1_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            }
                          }
                          else{
                            h_MJ1_um->Fill(jetpt, pt, dR, w_eec);// one is pythia
                            h3Jet_deltaR_MJ1_um->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ1_um->Fill(jetpt, dR, w_eec_3D);
                          }
                    }
                    else {
                      
                         h_MJ2->Fill(jetpt, pt, dR, w_eec);//both are pythia

                          h3Jet_deltaR_MJ2->Fill(jetpt, dR, w_eec_3D);
                          h3Jet_deltaR_MJ2->Fill(jetpt, dR, w_eec_3D);

                          if(tag=="matchedJet")
                          {
                            h_MJ2_m->Fill(jetpt, pt, dR, w_eec);// both is pythia
                            h_MJ2_tru_m->Fill(jetpt, pt, dR, w_eec_tru);// both is pythia

                            h3Jet_deltaR_MJ2_m->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ2_m->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MJ2_tru_m->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MJ2_tru_m->Fill(pt, dR, w_eec_tru_3D);

                            if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_MJ2_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MJ2_minpt->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MJ2_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MJ2_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            }
                          }
                         else{
                          h_MJ2_um->Fill(jetpt, pt, dR, w_eec);// both is pythia
                          h3Jet_deltaR_MJ2_um->Fill(jetpt, dR, w_eec_3D);
                          h3Jet_deltaR_MJ2_um->Fill(jetpt, dR, w_eec_3D);
                         }
                  }
                 }
                if(type == "sameMB")
                {
                  if(fCout)std::cout<<"In type: sameMB"<<endl;

                      h_MB1->Fill(jetpt, pt, dR, w_eec);
                      h3Jet_deltaR_MB1->Fill(jetpt, dR, w_eec_3D);
                      h3Jet_deltaR_MB1->Fill(jetpt, dR, w_eec_3D);
                     
                      if(tag=="matchedJet"){
                         h_MB1_m->Fill(jetpt, pt, dR, w_eec);
                         h_MB1_tru_m->Fill(jetpt, pt, dR, w_eec_tru);

                         h3Jet_deltaR_MB1_m->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_MB1_m->Fill(jetpt, dR, w_eec_3D);

                         h3Jet_deltaR_MB1_tru_m->Fill(pt, dR, w_eec_tru_3D);
                         h3Jet_deltaR_MB1_tru_m->Fill(pt, dR, w_eec_tru_3D);

                          if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_MB1_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MB1_minpt->Fill(jetpt, dR, w_eec_3D);

                            h3Jet_deltaR_MB1_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MB1_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                          }
                       }
                       else{
                            h_MB1_um->Fill(jetpt, pt, dR, w_eec);
                            h3Jet_deltaR_MB1_um->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MB1_um->Fill(jetpt, dR, w_eec_3D);
                       }
                }   
            }
        }
    }
    //else we are combining two different backgrounds or a jet + background
    else{
        int mult = particles.size();
        int mult2 = particles2.size();
        for (int i = 0; i < mult; i++){
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = 0; j < mult2; j++){
                 if(particles2.at(j).pt()<fMinENCtrackPt) continue;
                  
                   float pt_i = particles.at(i).pt();
                   float pt_j =  particles2.at(j).pt();
                    
                    w_eec = (1./(jetpt*jetpt))*(2*pt_i * pt_j);
                    w_eec_tru = (1./(pt*pt))*(2*pt_i * pt_j);

                    w_eec_3D =  (1.0 /(jetpt*jetpt))*(pt_i * pt_j);
                    w_eec_tru_3D =  (1.0 / (pt * pt))* (pt_i * pt_j);
                
                float dR = delR(particles.at(i), particles2.at(j));
                
                int i1 = particles.at(i).user_index();
                int i2 = particles2.at(j).user_index();
                
                if(type == "jetMB")
                {
                  if(fCout) cout<<"In jetMB"<<endl;
                    h_JMB->Fill(jetpt, pt, dR, w_eec);
                   
                    h3Jet_deltaR_JMB->Fill(jetpt, dR, w_eec_3D);
                    h3Jet_deltaR_JMB->Fill(jetpt, dR, w_eec_3D);

                    if(tag=="matchedJet"){
                      h_JMB_tru_m->Fill(jetpt, pt, dR, w_eec_tru);
                      h3Jet_deltaR_JMB_m->Fill(jetpt, dR, w_eec_3D);
                      h3Jet_deltaR_JMB_m->Fill(jetpt, dR, w_eec_3D);
                      h3Jet_deltaR_JMB_tru_m->Fill(pt, dR, w_eec_tru_3D);
                      h3Jet_deltaR_JMB_tru_m->Fill(pt, dR, w_eec_tru_3D); 

                      if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                        h3Jet_deltaR_JMB_minpt->Fill(jetpt, dR, w_eec_3D);
                        h3Jet_deltaR_JMB_minpt->Fill(jetpt, dR, w_eec_3D);
                      
                        h3Jet_deltaR_JMB_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                        h3Jet_deltaR_JMB_tru_minpt->Fill(pt, dR, w_eec_tru_3D); 
                      }
                    }
                    else{
                      h_JMB_um->Fill(jetpt, pt, dR, w_eec);
                      h3Jet_deltaR_JMB_um->Fill(jetpt, dR, w_eec_3D);
                      h3Jet_deltaR_JMB_um->Fill(jetpt, dR, w_eec_3D);
                    }

                  if(((i1 == bkgIndex1) && (i2 == bkgIndex2)) || ((i1 == bkgIndex2) && (i2 == bkgIndex1)))
                  {

                    if(fCout) cout<<"In BMB"<<endl;
                         h_BMB->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_BMB->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_BMB->Fill(jetpt, dR, w_eec_3D); 

                         if(tag=="matchedJet"){
                         h_BMB_tru_m->Fill(jetpt, pt, dR, w_eec_tru);
                         h_BMB_m->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_BMB_m->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_BMB_m->Fill(jetpt, dR, w_eec_3D); 

                         h3Jet_deltaR_BMB_tru_m->Fill(pt, dR, w_eec_tru_3D);
                         h3Jet_deltaR_BMB_tru_m->Fill(pt, dR, w_eec_tru_3D); 

                          if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_BMB_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_BMB_minpt->Fill(jetpt, dR, w_eec_3D); 

                            h3Jet_deltaR_BMB_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_BMB_tru_minpt->Fill(pt, dR, w_eec_tru_3D); 
                          }
                         }
                        else{
                         h_BMB_um->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_BMB_um->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_BMB_um->Fill(jetpt, dR, w_eec_3D); 

                        }
                  }
                   
                  else
                  {     
                    if(fCout) cout<<"In SMB"<<endl;
                         h_SMB->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_SMB->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_SMB->Fill(jetpt, dR, w_eec_3D); 

                         if(tag=="matchedJet"){
                         h_SMB_tru_m->Fill(jetpt, pt, dR, w_eec_tru);
                         h_SMB_m->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_SMB_m->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_SMB_m->Fill(jetpt, dR, w_eec_3D); 

                         h3Jet_deltaR_SMB_tru_m->Fill(pt, dR, w_eec_tru_3D);
                         h3Jet_deltaR_SMB_tru_m->Fill(pt, dR, w_eec_tru_3D); 

                          if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_SMB_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_SMB_minpt->Fill(jetpt, dR, w_eec_3D); 

                            h3Jet_deltaR_SMB_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_SMB_tru_minpt->Fill(pt, dR, w_eec_tru_3D); 
                          }
                         }
                        else{
                         h_SMB_um->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_SMB_um->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_SMB_um->Fill(jetpt, dR, w_eec_3D); 

                        }
                    }
                }
                if(type == "diffMB")
                {
                  if(fCout) cout<<"In diffMB"<<endl;
                  if(fCout) cout<<"distance "<<dR<<endl;
                    if ((i1 ==  bkgIndex1) && (i2 ==  bkgIndex2))
                    {
                    
                         h_MB1MB2->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_MB1MB2->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_MB1MB2->Fill(jetpt, dR, w_eec_3D); 

                         if(tag=="matchedJet"){
                         h_MB1MB2_tru_m->Fill(jetpt, pt, dR, w_eec_tru);
                         h_MB1MB2_m->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_MB1MB2_m->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_MB1MB2_m->Fill(jetpt, dR, w_eec_3D); 

                         h3Jet_deltaR_MB1MB2_tru_m->Fill(pt, dR, w_eec_tru_3D);
                         h3Jet_deltaR_MB1MB2_tru_m->Fill(pt, dR, w_eec_tru_3D); 

                          if(pt>fEmbTuthJetPtMin && ifMinPtHist){
                            h3Jet_deltaR_MB1MB2_minpt->Fill(jetpt, dR, w_eec_3D);
                            h3Jet_deltaR_MB1MB2_minpt->Fill(jetpt, dR, w_eec_3D); 

                            h3Jet_deltaR_MB1MB2_tru_minpt->Fill(pt, dR, w_eec_tru_3D);
                            h3Jet_deltaR_MB1MB2_tru_minpt->Fill(pt, dR, w_eec_tru_3D); 
                          }
                         }
                        else
                        {
                         h_MB1MB2_um->Fill(jetpt, pt, dR, w_eec); //background + MB

                         h3Jet_deltaR_MB1MB2_um->Fill(jetpt, dR, w_eec_3D);
                         h3Jet_deltaR_MB1MB2_um->Fill(jetpt, dR, w_eec_3D); 

                        }
                    }
                 }     
              }   
            }
        }
    if(fCout) cout<<"Finished with FillEmbJetsEEC"<<endl;
}

// //________________________________________________________________________
// //bkgindex = -1 for jet background particle, -2 for first thermal cone, -3 for second thermal cone, -4 for third thermal cone
// //pt is true (gen) level jet pT, jetpt is the pT of the embedded subtracted jet
void AliAnalysisTaskJetsEECpbpb::FillEmbJetsE3C(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2,std::vector<fastjet::PseudoJet> particles3, float jetpt, float pt, std::string typeSame, std::string type, bool ifMatchedJet)
{
    if(fCout){cout<<particles.size()<<" and "<<particles2.size()<<"  and "<<particles3.size()<<endl;}
    int mult = particles.size();
    int mult2 = particles2.size();
    int mult3 = particles3.size();
    
    float w = 0;
    float w_tru = 0;
    float w_iij = 0;
    float w_tru_iij = 0;
    float w_e3c_3D = 0;
    float w_e3c_tru_3D = 0;
    float w_iij_3D = 0;
    float w_iij_tru_3D = 0;
    
    float w_ijs = 0;
    float w_ijs_tru = 0;
    
    float w_ijs_3D = 0;
    float w_ijs_tru_3D = 0;
    
    float w_twosame = 0;
    float w_twosame_tru = 0;
    
    float w_twosame_3D = 0;
    float w_twosame_tru_3D = 0;
    
    float w_ijk = 0;
    float w_ijk_tru = 0;
    
    float w_ijk_3D = 0;
    float w_ijk_tru_3D = 0;
    
    
    if(fCout){cout<<"!!!!!!!!!!after setting: "<<mult<<" and "<<mult2<<" and "<<mult3<<endl;}
    if(fCout){cout<<"!!!!!!!!!!Setting are: typeSame "<<typeSame<<" and type "<<type<<endl;}
    
    if(typeSame == "all"){
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt() < fMinENCtrackPt ) continue;
            for (int j = i+1; j < mult; j++)
            {
                if(particles.at(j).pt()<fMinENCtrackPt) continue;
                
                if(fCout){cout<<"in this loop"<<endl;}
                
                
                w = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt());
                w_tru = (1./(pt*pt*pt))*(3*particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt());
                
                w_iij = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt());
                w_tru_iij = (1./(pt*pt*pt))*(3*particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt());
                
                
                w_e3c_3D =  (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt());
                w_e3c_tru_3D =  (1./(pt*pt*pt))*(particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt());
                
                w_iij_3D =  (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt());
                w_iij_tru_3D =  (1./(pt*pt*pt))*(particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt());
                
                
                
                float dR = delR(particles.at(i), particles.at(j));
                int i1 = particles.at(i).user_index();
                int i2 = particles.at(j).user_index();
                
                if(fCout){cout<<"weights calculated"<<endl;}
                
                if(type =="sameJet")
                {
                  if(fCout){cout<<"filling base histograms"<<endl;}
                    h_MJ_e3c->Fill(jetpt, pt, dR, w);h_MJ_e3c->Fill(jetpt, pt, dR, w_iij);

                  if(fCout){cout<<"filling base 3D histograms"<<endl;}

                    h3Jet_deltaR_MJ_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c->Fill(jetpt, dR, w_iij_3D);
                    h3Jet_deltaR_MJ_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c->Fill(jetpt, dR, w_iij_3D);
                    h3Jet_deltaR_MJ_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c->Fill(jetpt, dR, w_iij_3D);
                    
                    
                    if(ifMatchedJet){
                        if(fCout){cout<<"filling matched case"<<endl;}
                    
                        hJet_deltaR_MJ_e3c_m->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ_e3c_m->Fill(jetpt, pt, dR, w_iij);
                        hJet_deltaR_MJ_e3c_tru_m->Fill(jetpt, pt, dR, w_tru);hJet_deltaR_MJ_e3c_tru_m->Fill(jetpt, pt, dR, w_tru_iij);
                        
                        h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, dR, w_iij_3D);
                        h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, dR, w_iij_3D);
                        h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, dR, w_iij_3D);
                        
                       }
                        else{
                            if(fCout){cout<<"filling unmatched case"<<endl;}
                            hJet_deltaR_MJ_e3c_um->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ_e3c_um->Fill(jetpt, pt, dR, w_iij);
                            
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                  
                       }
                        
                        if ((i1 == -1) && (i2 == -1)) {
                          if(fCout){cout<<"all bkg"<<endl;}
                            h_MJ0_e3c->Fill(jetpt, pt, dR, w);h_MJ0_e3c->Fill(jetpt, pt, dR, w_iij);
                            
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c->Fill(jetpt, dR, w_iij_3D);
                            
                            if(ifMatchedJet){
                                
                                hJet_deltaR_MJ0_e3c_m->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ0_e3c_m->Fill(jetpt, pt, dR, w_iij);
                                hJet_deltaR_MJ0_e3c_tru_m->Fill(jetpt, pt, dR, w_tru);hJet_deltaR_MJ0_e3c_tru_m->Fill(jetpt, pt, dR, w_tru_iij);
                                
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                
                                
                            }
                            else{
                                
                                hJet_deltaR_MJ0_e3c_um->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ0_e3c_um->Fill(jetpt, pt, dR, w_iij);
                                
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                
                                
                            }
                        }//fake fake fake
                        else if (i2 == -1) {
                          if(fCout){cout<<"j bkg"<<endl;}
                            h_MJ1_e3c->Fill(jetpt, pt, dR, w);h_MJ2_e3c->Fill(jetpt, pt, dR, w_iij);

                            h3Jet_deltaR_MJ1_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ1_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ1_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c->Fill(jetpt, dR, w_iij_3D);
                            
                            if(ifMatchedJet){
                                
                                hJet_deltaR_MJ1_e3c_m->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ2_e3c_m->Fill(jetpt, pt, dR, w_iij);
                                hJet_deltaR_MJ1_e3c_tru_m->Fill(jetpt, pt, dR, w_tru);hJet_deltaR_MJ2_e3c_tru_m->Fill(jetpt, pt, dR, w_tru_iij);
                                
                                h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                
                                h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                
                            }
                            else{
                                
                                hJet_deltaR_MJ1_e3c_um->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ2_e3c_um->Fill(jetpt, pt, dR, w_iij);
                                
                                h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, dR, w_iij_3D);
                    
                            }
                        }//real fake fake for ijj & real real fake for iij
                        else if (i1 == -1) {
                          if(fCout){cout<<"i bkg, realrealfake"<<endl;}
                            h_MJ2_e3c->Fill(jetpt, pt, dR, w);h_MJ1_e3c->Fill(jetpt, pt, dR, w_iij);
                          
                            h3Jet_deltaR_MJ2_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ1_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ2_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ1_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ2_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ1_e3c->Fill(jetpt, dR, w_iij_3D);
                            
                            
                            if(ifMatchedJet){
                                
                                hJet_deltaR_MJ2_e3c_m->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ1_e3c_m->Fill(jetpt, pt, dR, w_iij);
                                hJet_deltaR_MJ2_e3c_tru_m->Fill(jetpt, pt, dR, w_tru);hJet_deltaR_MJ1_e3c_tru_m->Fill(jetpt, pt, dR, w_tru_iij);
                                
                                h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                
                            }
                            else{
                                
                                hJet_deltaR_MJ2_e3c_um->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ1_e3c_um->Fill(jetpt, pt, dR, w_iij);
                                
                                h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                
                                
                            }
                            
                        }//real real fake
                        else {
                            h_MJ3_e3c->Fill(jetpt, pt, dR, w);h_MJ3_e3c->Fill(jetpt, pt, dR, w_iij);
                            
                            h3Jet_deltaR_MJ3_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ3_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c->Fill(jetpt, dR, w_iij_3D);
                            h3Jet_deltaR_MJ3_e3c->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c->Fill(jetpt, dR, w_iij_3D);
                            
                            
                            if(ifMatchedJet){
                                hJet_deltaR_MJ3_e3c_m->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ3_e3c_m->Fill(jetpt, pt, dR, w_iij);
                                hJet_deltaR_MJ3_e3c_tru_m->Fill(jetpt, pt, dR, w_tru);hJet_deltaR_MJ3_e3c_tru_m->Fill(jetpt, pt, dR, w_tru_iij);
                                
                                h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, dR, w_iij_3D);
                                
                                h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, dR, w_iij_tru_3D);
                                
                            }
                            else{
                                hJet_deltaR_MJ3_e3c_um->Fill(jetpt, pt, dR, w);hJet_deltaR_MJ3_e3c_um->Fill(jetpt, pt, dR, w_iij);
                                
                                h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, dR, w_e3c_3D);h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, dR, w_iij_3D);
                                
                            }
                            
                        }//both are pythia

                        if(fCout){cout<<"finished with same jet loop"<<endl;}
                    }
                                    
                if(type=="sameMB"){
                    
                    if(fCout){cout<<"fill same MB"<<endl;}
                    h_MB1MB1MB1->Fill(jetpt, pt, dR, w);
                    h_MB1MB1MB1->Fill(jetpt, pt, dR, w_iij);
                    
                    h3_MB1MB1MB1->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1->Fill(jetpt, dR, w_iij_3D);
                    h3_MB1MB1MB1->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1->Fill(jetpt, dR, w_iij_3D);
                    h3_MB1MB1MB1->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1->Fill(jetpt, dR, w_iij_3D);
                    
                    
                    if(ifMatchedJet){
                        
                        h_MB1MB1MB1_m->Fill(jetpt, pt, dR, w);
                        h_MB1MB1MB1_m->Fill(jetpt, pt, dR, w_iij);
                        
                        h_MB1MB1MB1_tru_m->Fill(jetpt, pt, dR, w_tru);
                        h_MB1MB1MB1_tru_m->Fill(jetpt, pt, dR, w_tru_iij);
                        
                        h3_MB1MB1MB1_m->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_m->Fill(jetpt, dR, w_iij_3D);
                        h3_MB1MB1MB1_m->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_m->Fill(jetpt, dR, w_iij_3D);
                        h3_MB1MB1MB1_m->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_m->Fill(jetpt, dR, w_iij_3D);
                        
                        h3_MB1MB1MB1_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3_MB1MB1MB1_tru_m->Fill(pt, dR, w_iij_tru_3D);
                        h3_MB1MB1MB1_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3_MB1MB1MB1_tru_m->Fill(pt, dR, w_iij_tru_3D);
                        h3_MB1MB1MB1_tru_m->Fill(pt, dR, w_e3c_tru_3D);h3_MB1MB1MB1_tru_m->Fill(pt, dR, w_iij_tru_3D);
                        
                        
                        if(fCout){cout<<"matched same MB filled"<<endl;}
                    }
                    else{
                        
                        
                        h_MB1MB1MB1_um->Fill(jetpt, pt, dR, w);
                        h_MB1MB1MB1_um->Fill(jetpt, pt, dR, w_iij);
                        
                        h3_MB1MB1MB1_um->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_um->Fill(jetpt, dR, w_iij_3D);
                        h3_MB1MB1MB1_um->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_um->Fill(jetpt, dR, w_iij_3D);
                        h3_MB1MB1MB1_um->Fill(jetpt, dR, w_e3c_3D);h3_MB1MB1MB1_um->Fill(jetpt, dR, w_iij_3D);
                        
                        
                        
                    }
                    
                    if(fCout){cout<<"same MB filled"<<endl;}
                }
                for (int s = j+1; s < mult; s++)
                {
                    
                    if(particles.at(s).pt()<fMinENCtrackPt) continue;
                    
                    if(fCout){cout<<"s=j+1 loop"<<endl;}
                    
                    
                    w_ijs = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt());
                    w_ijs_tru = (1./(pt*pt*pt))*(6*particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt());
                    
                    w_ijs_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt());
                    w_ijs_tru_3D = (1./(pt*pt*pt))*(particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt());
                    
                    
                    float dR_ij = delR(particles.at(i), particles.at(j));
                    float dR_js = delR(particles.at(j), particles.at(s));
                    float dR_is = delR(particles.at(s), particles.at(i));
                    
                    float R_L = -1;
                    if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
                    else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
                    else{R_L = dR_is;}
                    
                    int i1 = particles.at(i).user_index();
                    int i2 = particles.at(j).user_index();
                    int i3 = particles.at(s).user_index();
                    
                    if(fCout){cout<<"ijs loop RL computed"<<endl;}
                    
                    if(type=="sameJet"){
                        
                        if(fCout){cout<<"fill same jet ijs"<<endl;}
                        
                        h_MJ_e3c->Fill(jetpt, pt, R_L, w_ijs);
                        
                        h3Jet_deltaR_MJ_e3c->Fill(jetpt, R_L, w_ijs_3D);
                        h3Jet_deltaR_MJ_e3c->Fill(jetpt, R_L, w_ijs_3D);
                        h3Jet_deltaR_MJ_e3c->Fill(jetpt, R_L, w_ijs_3D);
                        h3Jet_deltaR_MJ_e3c->Fill(jetpt, R_L, w_ijs_3D);
                        h3Jet_deltaR_MJ_e3c->Fill(jetpt, R_L, w_ijs_3D);
                        h3Jet_deltaR_MJ_e3c->Fill(jetpt, R_L, w_ijs_3D);
                        
                        if(ifMatchedJet){
                            
                            hJet_deltaR_MJ_e3c_m->Fill(jetpt, pt, R_L, w_ijs);
                            hJet_deltaR_MJ_e3c_tru_m->Fill(jetpt, pt, R_L, w_ijs_tru);
                            
                            h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                            
                            h3Jet_deltaR_MJ_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                            h3Jet_deltaR_MJ_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                            h3Jet_deltaR_MJ_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                            h3Jet_deltaR_MJ_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                            h3Jet_deltaR_MJ_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                            h3Jet_deltaR_MJ_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                            
                            
                        }
                        else{
                            
                            hJet_deltaR_MJ_e3c_um->Fill(jetpt, pt, R_L, w_ijs);
                            
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                            h3Jet_deltaR_MJ_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                            
                        }
                        
                        if ((i1 == -1) && (i2 == -1) && (i3 == -1)){
                            h_MJ0_e3c->Fill(jetpt, pt, R_L, w_ijs);
                            
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                            h3Jet_deltaR_MJ0_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                            
                            if(ifMatchedJet){
                                
                                hJet_deltaR_MJ0_e3c_m->Fill(jetpt, pt, R_L, w_ijs);
                                hJet_deltaR_MJ0_e3c_tru_m->Fill(jetpt, pt, R_L, w_ijs_tru);
                                
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake fake
                                
                                
                            }
                            else{
                                
                                hJet_deltaR_MJ0_e3c_um->Fill(jetpt, pt, R_L, w_ijs);
                                
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                h3Jet_deltaR_MJ0_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake fake
                                
                            }
                        }
                            else if (((i1 == -1) && (i2 == -1) && (i3 != 0)) ||((i1 == -1) && (i2 != 0) && (i3 == -1)) ||((i1 != 0) && (i2 == -1) && (i3 == -1)))
                            {
                                
                                h_MJ1_e3c->Fill(jetpt, pt, R_L, w_ijs);
                                
                                h3Jet_deltaR_MJ1_e3c->Fill(jetpt, R_L, w_ijs_3D);//fake fake REAL
                                h3Jet_deltaR_MJ1_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ1_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ1_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ1_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ1_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                
                                if(ifMatchedJet){
                                    
                                    
                                    hJet_deltaR_MJ1_e3c_m->Fill(jetpt, pt, R_L, w_ijs);
                                    hJet_deltaR_MJ1_e3c_tru_m->Fill(jetpt, pt, R_L, w_ijs_tru);
                                    
                                    h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, R_L, w_ijs_3D);//fake fake REAL
                                    h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    
                                    h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);//fake fake REAL
                                    h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ1_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    
                                    
                                }
                                else{
                                    
                                    hJet_deltaR_MJ1_e3c_um->Fill(jetpt, pt, R_L, w_ijs);
                                    
                                    h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, R_L, w_ijs_3D);//fake fake REAL
                                    h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ1_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    
                                    
                                }
                            }
                            else if (((i1 == -1) && (i2 != 0) && (i3 != 0)) ||((i1 != 0) && (i2 != 0) && (i3 == -1)) ||((i1 != 0) && (i2 == -1) && (i3 != 0)))
                            {
                                h_MJ2_e3c->Fill(jetpt, pt, R_L, w_ijs);
                                
                                h3Jet_deltaR_MJ2_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ2_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ2_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ2_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ2_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ2_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                
                                if(ifMatchedJet){
                                    
                                    hJet_deltaR_MJ2_e3c_m->Fill(jetpt, pt, R_L, w_ijs);
                                    hJet_deltaR_MJ2_e3c_tru_m->Fill(jetpt, pt, R_L, w_ijs_tru);
                                    
                                    h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    
                                    h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ2_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    
                                    
                                }
                                else{
                                    
                                    hJet_deltaR_MJ2_e3c_um->Fill(jetpt, pt, R_L, w_ijs);
                                    
                                    h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ2_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    
                                }
                                
                            }
                            else{
                                h_MJ3_e3c->Fill(jetpt, pt, R_L, w_ijs);//real real real
                                
                                
                                h3Jet_deltaR_MJ3_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ3_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ3_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ3_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ3_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                h3Jet_deltaR_MJ3_e3c->Fill(jetpt, R_L, w_ijs_3D);
                                
                                if(ifMatchedJet){
                                    
                                    
                                    hJet_deltaR_MJ3_e3c_m->Fill(jetpt, pt, R_L, w_ijs);//real real real
                                    hJet_deltaR_MJ3_e3c_tru_m->Fill(jetpt, pt, R_L, w_ijs_tru);
                                    
                                    h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_m->Fill(jetpt, R_L, w_ijs_3D);
                                    
                                    h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    h3Jet_deltaR_MJ3_e3c_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                    
                                    
                                }
                                else{
                                    
                                    hJet_deltaR_MJ3_e3c_um->Fill(jetpt, pt, R_L, w_ijs);//real real real
                                    
                                    h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    h3Jet_deltaR_MJ3_e3c_um->Fill(jetpt, R_L, w_ijs_3D);
                                    
                                }
                            }
                        }
                        
                        if(type=="sameMB"){
                            if(fCout){cout<<"same MB ijs loop"<<endl;}
                            h_MB1MB1MB1->Fill(jetpt, pt, R_L, w_ijs);
                            
                            h3_MB1MB1MB1->Fill(jetpt, R_L, w_ijs_3D);
                            h3_MB1MB1MB1->Fill(jetpt, R_L, w_ijs_3D);
                            h3_MB1MB1MB1->Fill(jetpt, R_L, w_ijs_3D);
                            h3_MB1MB1MB1->Fill(jetpt, R_L, w_ijs_3D);
                            h3_MB1MB1MB1->Fill(jetpt, R_L, w_ijs_3D);
                            h3_MB1MB1MB1->Fill(jetpt, R_L, w_ijs_3D);
                            
                            if(ifMatchedJet){
                                
                                h_MB1MB1MB1_m->Fill(jetpt, pt, R_L, w_ijs);
                                h_MB1MB1MB1_tru_m->Fill(jetpt, pt, R_L, w_ijs_tru);
                                
                                h3_MB1MB1MB1_m->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_m->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_m->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_m->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_m->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_m->Fill(jetpt, R_L, w_ijs_3D);
                                
                                h3_MB1MB1MB1_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                h3_MB1MB1MB1_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                h3_MB1MB1MB1_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                h3_MB1MB1MB1_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                h3_MB1MB1MB1_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                h3_MB1MB1MB1_tru_m->Fill(pt, R_L, w_ijs_tru_3D);
                                
                            }
                            else{
                                
                                h_MB1MB1MB1_um->Fill(jetpt, pt, R_L, w_ijs);
                                
                                h3_MB1MB1MB1_um->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_um->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_um->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_um->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_um->Fill(jetpt, R_L, w_ijs_3D);
                                h3_MB1MB1MB1_um->Fill(jetpt, R_L, w_ijs_3D);
                                
                                
                            }
                        }
                    }
                }
            }
        }
    
    else if (typeSame == "two"){
        //particles2 and particles3 is the same list
        
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            int i1 = particles.at(i).user_index();
            
            for (int  j = 0; j < mult2; j++)
            {
                
                if(particles2.at(j).pt()<fMinENCtrackPt) continue;
                
                int i2 = particles2.at(j).user_index();
                
                w_twosame = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt());
                w_twosame_tru = (1./(pt*pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt());
                
                w_twosame_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt());
                w_twosame_tru_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt());
                
                float R_L = delR(particles.at(i), particles2.at(j));
                
                if(fCout){cout<<" about to fill histograms in TWO SAME"<<endl;}
                if(i1==-2 && (i2 >= 0 || i2 == -1)) {
                    
                    if(fCout){cout<<" filled jjmb TWO SAME"<<endl;}
                    
                    h_JJMB->Fill(jetpt, pt, R_L,w_twosame);
                    
                    h3_JJMB->Fill(jetpt, pt, R_L,w_twosame_3D);
                    h3_JJMB->Fill(jetpt, pt, R_L,w_twosame_3D);
                    h3_JJMB->Fill(jetpt, pt, R_L,w_twosame_3D);
                    
                    if(ifMatchedJet){
                        h_JJMB_m->Fill(jetpt, pt, R_L,w_twosame);
                        h_JJMB_tru_m->Fill(jetpt, pt, R_L,w_twosame_tru);
                        
                        h3_JJMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_JJMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_JJMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        
                        h3_JJMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_JJMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_JJMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                    }
                    else{
                        
                        h_JJMB_um->Fill(jetpt, pt, R_L,w_twosame);
                        
                        h3_JJMB_um->Fill(jetpt, pt, R_L,w_twosame_3D);
                        h3_JJMB_um->Fill(jetpt, pt, R_L,w_twosame_3D);
                        h3_JJMB_um->Fill(jetpt, pt, R_L,w_twosame_3D);
                        
                        
                    }
                }
                if(i1==-2 && i2 == -1 && i2 == -1){
                    if(fCout){cout<<" filled bbmb TWO SAME"<<endl;}
                    
                    h_BBMB->Fill(jetpt, pt, R_L,w_twosame);
                    
                    h3_BBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_BBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_BBMB->Fill(jetpt, R_L,w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        h_BBMB_m->Fill(jetpt, pt, R_L,w_twosame);
                        h_BBMB_tru_m->Fill(jetpt, pt, R_L,w_twosame_tru);
                        
                        h3_BBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        
                        h3_BBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_BBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_BBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        
                    }
                    else{
                        
                        h_BBMB_um->Fill(jetpt, pt, R_L,w_twosame);
                        
                        h3_BBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        
                    }
                }
                if(i1==-2 && i2 >= 0 && i2 >= 0){
                    
                    if(fCout){cout<<" filled ssmb TWO SAME"<<endl;}
                    
                    h_SSMB->Fill(jetpt, pt, R_L,w_twosame);
                    
                    h3_SSMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_SSMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_SSMB->Fill(jetpt, R_L,w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        h_SSMB_m->Fill(jetpt, pt, R_L,w_twosame);
                        h_SSMB_tru_m->Fill(jetpt, pt, R_L,w_twosame_tru);
                        
                        h3_SSMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SSMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SSMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        
                        h3_SSMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_SSMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_SSMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        
                        
                    }
                    else{
                        
                        h_SSMB_um->Fill(jetpt, pt, R_L,w_twosame);
                        
                        h3_SSMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SSMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SSMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        
                    }
                }
                
                if((i1 == -1 || i1 >= 0) && (i2 == -2)){
                    
                    if(fCout){cout<<" filled jmbmb TWO SAME"<<endl;}
                    
                    h_JMBMB->Fill(jetpt, pt, R_L,w_twosame);
                    
                    h3_JMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_JMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_JMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        h_JMBMB_m->Fill(jetpt, pt, R_L,w_twosame);
                        h_JMBMB_tru_m->Fill(jetpt, pt, R_L,w_twosame_tru);
                        
                        h3_JMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_JMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_JMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        
                        h3_JMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_JMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_JMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        
                    }
                    else{
                        
                        h_JMBMB_um->Fill(jetpt, pt, R_L,w_twosame);
                        
                        h3_JMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_JMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_JMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        
                    }
                }
                if((i1== -1)  && (i2 == -2) && (i2==-2)){
                    
                    if(fCout){cout<<" filled bmbmb TWO SAME"<<endl;}
                    
                    h_BMBMB->Fill(jetpt, pt, R_L,w_twosame);
                    
                    h3_BMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_BMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_BMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        h_BMBMB_m->Fill(jetpt, pt, R_L,w_twosame);
                        h_BMBMB_tru_m->Fill(jetpt, pt, R_L,w_twosame_tru);
                        
                        h3_BMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        
                        h3_BMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_BMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_BMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        
                    }
                    else{
                        
                        h_BMBMB_um->Fill(jetpt, pt, R_L,w_twosame);
                        
                        h3_BMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_BMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        
                        
                    }
                }
                if((i1 >= 0) && (i2 == -2) && (i2==-2)){
                    
                    if(fCout){cout<<" filled smbmb TWO SAME"<<endl;}
                    
                    h_SMBMB->Fill(jetpt, pt, R_L,w_twosame);
                    
                    h3_SMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_SMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    h3_SMBMB->Fill(jetpt, R_L,w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        h_SMBMB_m->Fill(jetpt, pt, R_L,w_twosame);
                        h_SMBMB_tru_m->Fill(jetpt, pt, R_L,w_twosame_tru);
                        
                        h3_SMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SMBMB_m->Fill(jetpt, R_L,w_twosame_3D);
                        
                        h3_SMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_SMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        h3_SMBMB_tru_m->Fill(pt, R_L,w_twosame_tru_3D);
                        
                    }
                    else{
                        
                        h_SMBMB_um->Fill(jetpt, pt, R_L,w_twosame);
                        
                        h3_SMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        h3_SMBMB_um->Fill(jetpt, R_L,w_twosame_3D);
                        
                    }
                }
                
                if(i1 == -3 && i2 == -2){
                    
                    if(fCout){cout<<" filled mb1mb1mb2 TWO SAME"<<endl;}
                    
                    h_MB1MB1MB2->Fill(jetpt, pt, R_L, w_twosame);
                    
                    h3_MB1MB1MB2->Fill(jetpt, R_L, w_twosame_3D);
                    h3_MB1MB1MB2->Fill(jetpt, R_L, w_twosame_3D);
                    h3_MB1MB1MB2->Fill(jetpt, R_L, w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        h_MB1MB1MB2_m->Fill(jetpt, pt, R_L, w_twosame);
                        h_MB1MB1MB2_tru_m->Fill(jetpt, pt, R_L, w_twosame_tru);
                        
                        h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_twosame_3D);
                        
                        h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_twosame_tru_3D);
                        h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_twosame_tru_3D);
                        h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_twosame_tru_3D);
                        
                    }
                    else{
                        
                        h_MB1MB1MB2_um->Fill(jetpt, pt, R_L, w_twosame);
                        
                        h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_twosame_3D);
                        
                    }
                }
                if(i1 == -2 && i2 == -3){
                    
                    if(fCout){cout<<" filled mb1mb2mb2 TWO SAME"<<endl;}
                    
                    h_MB1MB2MB2->Fill(jetpt, pt, R_L, w_twosame);
                    
                    h3_MB1MB2MB2->Fill(jetpt, R_L, w_twosame_3D);
                    h3_MB1MB2MB2->Fill(jetpt, R_L, w_twosame_3D);
                    h3_MB1MB2MB2->Fill(jetpt, R_L, w_twosame_3D);
                    
                    if(ifMatchedJet){
                        
                        
                        h_MB1MB2MB2_m->Fill(jetpt, pt, R_L, w_twosame);
                        h_MB1MB2MB2_tru_m->Fill(jetpt, pt, R_L, w_twosame_tru);
                        
                        h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_twosame_3D);
                        
                        h3_MB1MB2MB2_tru_m->Fill(jetpt, R_L, w_twosame_tru_3D);
                        h3_MB1MB2MB2_tru_m->Fill(jetpt, R_L, w_twosame_tru_3D);
                        h3_MB1MB2MB2_tru_m->Fill(jetpt, R_L, w_twosame_tru_3D);
                        
                    }
                    else{
                        
                        h_MB1MB2MB2_um->Fill(jetpt, pt, R_L, w_twosame);
                        
                        h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_twosame_3D);
                        h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_twosame_3D);
                        
                    }
                }
                
                if(fCout){cout<<" filled j loop in TWO SAME"<<endl;}
                
                for (int k = j+1; k < mult2; k++)
                {
                    if(particles2.at(k).pt()<fMinENCtrackPt) continue;
                    int i3 = particles2.at(k).user_index();
                    
                    
                    w_ijk = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt());
                    w_ijk_tru = (1./(pt*pt*pt))*(6*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt());
                    
                    w_ijk_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt());
                    w_ijk_tru_3D = (1./(pt*pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt());
                    
                    
                    
                    float dR_ij = delR(particles.at(i), particles2.at(j));
                    float dR_js = delR(particles2.at(j), particles2.at(k));
                    float dR_is = delR(particles2.at(k), particles.at(i));
                    
                    float R_L = -1;
                    if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
                    else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
                    else{R_L = dR_is;}
                    
                    
                    if(i1==-2 && (i2 >= 0 || i2 == -1) && (i3 >= 0 || i3 == -1)) {
                        
                        if(fCout){cout<<" filled jjmb two same "<<endl;}
                        h_JJMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_JJMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JJMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JJMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JJMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JJMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JJMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_JJMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_JJMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_JJMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_JJMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_JJMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_JJMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_JJMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_JJMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_JJMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_JJMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_JJMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JJMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                    if(i1==-2 && i2 == -1 && i3 == -1){
                        if(fCout){cout<<" filled bbmb two same "<<endl;}
                        
                        h_BBMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_BBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BBMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_BBMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_BBMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_BBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_BBMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_BBMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_BBMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_BBMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_BBMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            h3_BBMB_tru_m->Fill(jetpt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_BBMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_BBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                        }
                    }
                    if(i1==-2 && i2 >= 0 && i3 >= 0){
                        
                        if(fCout){cout<<" filled ssmb two same "<<endl;}
                        
                        h_SSMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_SSMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SSMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SSMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SSMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SSMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SSMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_SSMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_SSMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_SSMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_SSMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SSMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SSMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SSMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SSMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SSMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_SSMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_SSMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SSMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                    if(i1==-2 && ((i2 == -1 && i3 >= 0)||(i2 >= 0 && i3 == -1))) {
                        
                        if(fCout){cout<<" filled sbmb two same "<<endl;}
                        h_SBMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_SBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SBMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_SBMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_SBMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_SBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_SBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_SBMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_SBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                        }
                    }
                    if((i1== -1 || i1 >= 0) && (i2 == -2) && (i3==-2)) {
                        
                        if(fCout){cout<<" filled jmbmb two same "<<endl;}
                        h_JMBMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_JMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_JMBMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_JMBMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_JMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_JMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_JMBMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_JMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                    if((i1== -1)  && (i2 == -2) && (i3== -2)){
                        
                        if(fCout){cout<<" filled bmbmb two same "<<endl;}
                        h_BMBMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_BMBMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_BMBMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            
                            h_BMBMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                            
                        }
                    }
                    if((i1 >= 0) && (i2 == -2) && (i3== -2)){
                        
                        if(fCout){cout<<" filled smbmb two same "<<endl;}
                        
                        h_SMBMB->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_SMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMBMB->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_SMBMB_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_SMBMB_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_SMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_SMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMBMB_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_SMBMB_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_SMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMBMB_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                        }
                    }
                    
                    if(i1 == -3 && i2 == -2 && i3 == -2){
                        
                        if(fCout){cout<<" filled mb1mb1mb2 two same "<<endl;}
                        
                        h_MB1MB1MB2->Fill(jetpt, pt, R_L, w_ijk);
                        
                        h3_MB1MB1MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB1MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB1MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB1MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB1MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB1MB2->Fill(jetpt, R_L, w_ijk_3D);
                        if(ifMatchedJet){
                            
                            h_MB1MB1MB2_m->Fill(jetpt, pt, R_L, w_ijk);
                            h_MB1MB1MB2_tru_m->Fill(jetpt, pt, R_L, w_ijk_tru);
                            
                            h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            
                            h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB1MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            
                            
                        }
                        else{
                            
                            h_MB1MB1MB2_um->Fill(jetpt, pt, R_L, w_ijk);
                            
                            h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB1MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            
                        }
                    }
                    if(i1 == -2 && i2 == -3 && i3 == -3){
                        if(fCout){cout<<" filled mb1mb2mb2 two same "<<endl;}
                        
                        h_MB1MB2MB2->Fill(jetpt, pt, R_L, w_ijk);
                        
                        h3_MB1MB2MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB2MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB2MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB2MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB2MB2->Fill(jetpt, R_L, w_ijk_3D);
                        h3_MB1MB2MB2->Fill(jetpt, R_L, w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_MB1MB2MB2_m->Fill(jetpt, pt, R_L, w_ijk);
                            h_MB1MB2MB2_tru_m->Fill(jetpt, pt, R_L, w_ijk_tru);
                            
                            h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_m->Fill(jetpt, R_L, w_ijk_3D);
                            
                            h3_MB1MB2MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB2MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB2MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB2MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB2MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            h3_MB1MB2MB2_tru_m->Fill(pt, R_L, w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_MB1MB2MB2_um->Fill(jetpt, pt, R_L, w_ijk);
                            
                            h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            h3_MB1MB2MB2_um->Fill(jetpt, R_L, w_ijk_3D);
                            
                        }
                        
                    }
                }
            }
        }
    }
    else{
        if(fCout){cout<<" all diff now "<<endl;}
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            int i1 = particles.at(i).user_index();
            
            for (int  j = 0; j < mult2; j++)
            {
                if(particles2.at(j).pt()<fMinENCtrackPt) continue;
                int i2 = particles2.at(j).user_index();
                
                for (int k = 0; k < mult3; k++)
                {
                    if(particles3.at(k).pt()<fMinENCtrackPt) continue;
                    int i3 = particles3.at(k).user_index();
                    
                    w_ijk = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt());
                    w_ijk_tru = (1./(pt*pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt());
                    
                    w_ijk_3D = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt());
                    w_ijk_tru_3D = (1./(pt*pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt());
                    
                    float dR_ij = delR(particles.at(i), particles2.at(j));
                    float dR_js = delR(particles2.at(j), particles3.at(k));
                    float dR_is = delR(particles.at(i), particles3.at(k));
                    
                    float R_L = -1;
                    if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
                    else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
                    else{R_L = dR_is;}
                    
                    if((i1 >= 0 || i1 == -1) && i2==-2 && i3 == -3){
                        if(fCout){cout<<" jmb1mb2 "<<endl;}
                        
                        h_JMB1MB2->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_JMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_JMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_JMB1MB2_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_JMB1MB2_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_JMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_JMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_JMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_JMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk_3D);
                            h3_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk_3D);
                            h3_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk_3D);
                            h3_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk_3D);
                            h3_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk_3D);
                            h3_JMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                    if(i1 == -1 && i2==-2 && i3 == -3){
                        
                        if(fCout){cout<<" jbmb1mb2 "<<endl;}
                        
                        h_BMB1MB2->Fill(jetpt, pt, R_L,w_ijk);
                        h_BMB1MB2_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
                        
                        h3_BMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_BMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            
                            h_BMB1MB2_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_BMB1MB2_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_BMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_BMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_BMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                        }
                        
                        else{
                            
                            h_BMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_BMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_BMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                    if(i1 >= 0 && i2==-2 && i3 == -3){
                        
                        if(fCout){cout<<" sbmb1mb2 "<<endl;}
                        
                        h_SMB1MB2->Fill(jetpt, pt, R_L,w_ijk);
                        
                        h3_SMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        h3_SMB1MB2->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_SMB1MB2_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_SMB1MB2_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_SMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_SMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_SMB1MB2_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_SMB1MB2_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_SMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_SMB1MB2_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                    
                    if(i1 == -2 && i2==-3 && i3 == -4){
                        
                        if(fCout){cout<<" mb1mb2mb3 "<<endl;}
                        
                        h_MB1MB2MB3->Fill(jetpt, pt, R_L,w_ijk);
                        h_MB1MB2MB3_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
                        
                        h3_MB1MB2MB3->Fill(jetpt, R_L,w_ijk_3D);
                        h3_MB1MB2MB3->Fill(jetpt, R_L,w_ijk_3D);
                        h3_MB1MB2MB3->Fill(jetpt, R_L,w_ijk_3D);
                        h3_MB1MB2MB3->Fill(jetpt, R_L,w_ijk_3D);
                        h3_MB1MB2MB3->Fill(jetpt, R_L,w_ijk_3D);
                        h3_MB1MB2MB3->Fill(jetpt, R_L,w_ijk_3D);
                        
                        if(ifMatchedJet){
                            
                            h_MB1MB2MB3_m->Fill(jetpt, pt, R_L,w_ijk);
                            h_MB1MB2MB3_tru_m->Fill(jetpt, pt, R_L,w_ijk_tru);
                            
                            h3_MB1MB2MB3_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_m->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_m->Fill(jetpt, R_L,w_ijk_3D);
                            
                            h3_MB1MB2MB3_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_MB1MB2MB3_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_MB1MB2MB3_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_MB1MB2MB3_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_MB1MB2MB3_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            h3_MB1MB2MB3_tru_m->Fill(pt, R_L,w_ijk_tru_3D);
                            
                        }
                        else{
                            
                            h_MB1MB2MB3_um->Fill(jetpt, pt, R_L,w_ijk);
                            
                            h3_MB1MB2MB3_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_um->Fill(jetpt, R_L,w_ijk_3D);
                            h3_MB1MB2MB3_um->Fill(jetpt, R_L,w_ijk_3D);
                            
                            
                        }
                    }
                }
            }
        }
    }
}
// //______________________________________________________________________
// //for embedded, det level and truth jets 
std::pair<std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>>  AliAnalysisTaskJetsEECpbpb::FindThermalConeEEC(AliEmcalJet *fJetEmb, int index1, int index2)
{
    std::vector<fastjet::PseudoJet> coneParticlesPlus;
    std::vector<fastjet::PseudoJet> coneParticlesMinus;

    Float_t jet_embpt = fJetEmb->Pt();
    Float_t jet_embphi = fJetEmb->Phi();
    Float_t jet_embeta = fJetEmb->Eta();
    Float_t jet_embArea = fJetEmb->AreaPt();

    Double_t dPhi1 = 999.;
    Double_t dPhi2 = 999.;
    Double_t dPhi = 999.;
    Double_t dEta = 999.;
    Double_t Axis1 = 999, Axis2 = 999;
    Double_t fEtaMC = 999;
    fastjet::PseudoJet PseudoTracks; // Creating a pseudojet object

    if (fDoRandCone){
      double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // Random Phi in the range of π/3 and 2π/3
      Axis1 = jet_embphi + OffsetRandom;
      Axis2 = jet_embphi - OffsetRandom;
      if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi();
      if(Axis2 < 0) Axis2 += TMath::TwoPi();
    }

    if (fDoPerpCone){
      Axis1 = ((jet_embphi + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? jet_embphi - ((3. / 2.) * TMath::Pi()) : jet_embphi + (TMath::Pi() / 2.);
      Axis2 = ((jet_embphi - (TMath::Pi() / 2.)) < 0) ? jet_embphi + ((3. / 2.) * TMath::Pi()) : jet_embphi - (TMath::Pi() / 2.);
    }

    AliParticleContainer *partCont = 0;
    TIter nextPartCont(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer *>(nextPartCont()))) {
      AliParticleIterableMomentumContainer itcont = partCont->accepted_momentum();
      for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
        AliVTrack * particle = static_cast<AliVTrack*>(it->second);
        AliAODTrack* trackReal = (AliAODTrack*)(particle);
        
            if (!trackReal) continue;
            if (TMath::Abs(trackReal->Eta()) > fEtaCutValue) continue;
            if (trackReal->Pt() < fMinENCtrackPt || trackReal->Pt() > fMaxPtTrack ) continue;

            Float_t track_phi = trackReal->Phi();
            Float_t track_eta = trackReal->Eta();

            Double_t dPhi1 = TMath::Abs(track_phi - Axis1);
            dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
            Double_t dPhi2 = TMath::Abs(track_phi - Axis2);
            dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
            Double_t dEta = jet_embeta - track_eta;

            if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) < fConeR)) {
                Float_t rotated_phi = jet_embphi + (track_phi - Axis1);
                if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
                if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();

                fastjet::PseudoJet pseudoTrack(trackReal->Pt() * TMath::Cos(rotated_phi), 
                                               trackReal->Pt() * TMath::Sin(rotated_phi),
                                               trackReal->Pt() * TMath::SinH(track_eta),
                                               trackReal->E());
                pseudoTrack.set_user_index(index1);
                coneParticlesPlus.push_back(pseudoTrack);
            }

            if ((TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) < fConeR)) {
                Float_t rotated_phi = jet_embphi + (track_phi - Axis2);
                if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
                if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();

                fastjet::PseudoJet pseudoTrack(trackReal->Pt() * TMath::Cos(rotated_phi), 
                                               trackReal->Pt() * TMath::Sin(rotated_phi),
                                               trackReal->Pt() * TMath::SinH(track_eta),
                                               trackReal->E());
                pseudoTrack.set_user_index(index2);
                coneParticlesMinus.push_back(pseudoTrack);
            }
        }
    }
    return std::make_pair(coneParticlesPlus, coneParticlesMinus);
}
// //______________________________________________________________________
// //outputs three cones 
std::tuple<std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>>
AliAnalysisTaskJetsEECpbpb::FindThermalConeE3C(AliEmcalJet *fJetEmb, int index1, int index2, int index3) 
{
    std::vector<fastjet::PseudoJet> coneParticles1, coneParticles2, coneParticles3; // Separate vectors for each cone

    // Jet kinematics
    Float_t jet_embphi = fJetEmb->Phi();
    Float_t jet_embeta = fJetEmb->Eta();

    // Anti-jet location (π - jet_embphi)
    Float_t anti_jet_phi = jet_embphi + TMath::Pi();
    if (anti_jet_phi > TMath::TwoPi()) anti_jet_phi -= TMath::TwoPi();

    Double_t Axis1_Perp, Axis2_Shifted, Axis3_Offset;
    Double_t dPhi1, dPhi2, dPhi3, dEta1, dEta2, dEta3;
    Double_t fEtaMC;
    fastjet::PseudoJet PseudoTracks;

    // **Cone 1: Perpendicular to jet axis (π/2 shift)**
    Axis1_Perp = jet_embphi + (TMath::Pi() / 2.);
    if (Axis1_Perp > TMath::TwoPi()) Axis1_Perp -= TMath::TwoPi();

    // **Cone 2: Shift φ by fDeltaAxisShift radians and reflect in η**
    Axis2_Shifted = jet_embphi + fDeltaAxisShift;
    if (Axis2_Shifted > TMath::TwoPi()) Axis2_Shifted -= TMath::TwoPi();
    Float_t eta_reflected = -jet_embeta; // Reflect in η for selection

    // **Cone 3: fDeltaAxisShift away from Cone 2, same η as the jet**
    Axis3_Offset = Axis2_Shifted + fDeltaAxisShift;
    if (Axis3_Offset > TMath::TwoPi()) Axis3_Offset -= TMath::TwoPi();

    AliParticleContainer *partCont = 0;
    TIter nextPartCont(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer *>(nextPartCont()))) {
        AliParticleIterableMomentumContainer itcont = partCont->accepted_momentum();
        for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
            AliVTrack *particle = static_cast<AliVTrack *>(it->second);
            AliAODTrack *trackReal = (AliAODTrack *)(particle);

            if (!trackReal) continue;

            fEtaMC = trackReal->Eta();
            if (TMath::Abs(fEtaMC) > fEtaCutValue) continue;
            if (trackReal->Pt() < fMinENCtrackPt || trackReal->Pt() > fMaxPtTrack) continue;

            Float_t track_phi = trackReal->Phi();
            Float_t track_eta = trackReal->Eta();

            // Compute distances to cones using original η-φ values
            dPhi1 = TMath::Abs(track_phi - Axis1_Perp);
            dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
            dEta1 = jet_embeta - track_eta;
            double distanceCone1 = TMath::Sqrt(dPhi1 * dPhi1 + dEta1 * dEta1);

            dPhi2 = TMath::Abs(track_phi - Axis2_Shifted);
            dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
            dEta2 = eta_reflected - track_eta;
            double distanceCone2 = TMath::Sqrt(dPhi2 * dPhi2 + dEta2 * dEta2);

            dPhi3 = TMath::Abs(track_phi - Axis3_Offset);
            dPhi3 = (dPhi3 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi3 : dPhi3;
            dEta3 = jet_embeta - track_eta;
            double distanceCone3 = TMath::Sqrt(dPhi3 * dPhi3 + dEta3 * dEta3);


            // Assign to correct cone
            if (distanceCone1 < fConeR) {
              Float_t rotated_phi = jet_embphi + (track_phi - Axis1_Perp);  // Rotation in φ
              if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
              if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
              fastjet::PseudoJet PseudoTracks(trackReal->Pt() * TMath::Cos(rotated_phi), 
                                          trackReal->Pt() * TMath::Sin(rotated_phi),
                                          trackReal->Pt() * TMath::SinH(track_eta), // Keeping η the same
                                          trackReal->E());
              PseudoTracks.set_user_index(index1);  
              coneParticles1.push_back(PseudoTracks);
            }
            if (distanceCone2 < fConeR) {
              Float_t rotated_phi = jet_embphi + (track_phi - Axis2_Shifted);  // Rotation in φ
              if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
              if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
              Float_t new_eta = -track_eta;
              fastjet::PseudoJet PseudoTracks(trackReal->Pt() * TMath::Cos(rotated_phi), 
                                          trackReal->Pt() * TMath::Sin(rotated_phi),
                                          trackReal->Pt() * TMath::SinH(new_eta), // Keeping η the same
                                          trackReal->E());
              PseudoTracks.set_user_index(index2); // Assign index to track            
              coneParticles2.push_back(PseudoTracks);
            }
            if (distanceCone3 < fConeR) {
              Float_t rotated_phi = jet_embphi + (track_phi - Axis3_Offset);  // Rotation in φ
              if (rotated_phi > TMath::Pi()) rotated_phi -= 2 * TMath::Pi();
              if (rotated_phi < -TMath::Pi()) rotated_phi += 2 * TMath::Pi();
              fastjet::PseudoJet PseudoTracks(trackReal->Pt() * TMath::Cos(rotated_phi), 
                                          trackReal->Pt() * TMath::Sin(rotated_phi),
                                          trackReal->Pt() * TMath::SinH(track_eta), // Keeping η the same
                                          trackReal->E());
              PseudoTracks.set_user_index(index3); // Assign index to track           
              coneParticles3.push_back(PseudoTracks);
            }
        }
    }

    return std::make_tuple(coneParticles1, coneParticles2, coneParticles3);
}

// //______________________________________________________________________   
//ENC computation-------------------------------------------------------
void AliAnalysisTaskJetsEECpbpb::ComputeENC(AliEmcalJet *fJet, float ptSub, AliJetContainer *fJetCont)
{
    //General ENC computation: Need to loop over jets, then loop within a single jet to compute EECs.
    //NOTE: Already in the jet loop defined in the Fill Histogram function. This puts you in the jet loop. The event loop stuff is being taken care of by AliAnalysis manager
    //fjet is my single jet. AliEmCalJet is pointing to the fJet object.
    
    std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
    fConstituents.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex = 0;
    //The line below gets constituent particles within fjet. C++ syntax[ for (auto elem : container)    // capture elements by value ]
    for (auto part: fJet->GetParticleConstituents())
    {
        PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet)); //leave this in for the same reason as above
        if (PseudoTracks.pt() < fMinENCtrackPt || PseudoTracks.pt() > fMaxPtTrack) continue; //remove tracks below cut for ENCs
        fConstituents.push_back(PseudoTracks);
        constituentIndex++;
    }
    
    if(fCout)cout<<"In Compute ENC"<<endl;
       
        //Looping over the jet
        // double jet_pt = fJet->Pt();
        float jet_pt = ptSub;
        jet_pt_hist->Fill(jet_pt);
        if(fCout)cout<<"Got subtracted jet pT"<<endl;
          //Looping over the det jet
     for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
        {
            for(int s=j+1; s<int(fConstituents.size()) ; s++)
            {
              // if(fCout)cout<<"in loop s"<<endl;

                float eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                float eee_jjs =((3*fConstituents[j].pt()*fConstituents[j].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));

                float eee_jss_3D =((fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                float eee_jjs_3D =((fConstituents[j].pt()*fConstituents[j].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));

                float deltaR = delR(fConstituents[j],fConstituents[s]);

                float ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                float ee_js_3D = (fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));

                if(fpaircut == 1)
                {
                    float j_eta = fConstituents[j].eta();
                    float s_eta = fConstituents[s].eta();
                    float del_js_eta = abs(j_eta-s_eta);
                    if (del_js_eta < 0.008) continue;
                    else
                    {

                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jss);
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jjs);

                        EEC_pt_hist->Fill(deltaR,jet_pt,ee_js);

                        OptUn_eec->Fill(jet_pt,deltaR,ee_js_3D);
                        OptUn_eec->Fill(jet_pt,deltaR,ee_js_3D);

                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jss_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jss_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jss_3D);

                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jjs_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jjs_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jjs_3D);


                    }
                }
                else
                {
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jss);
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jjs);

                        EEC_pt_hist->Fill(deltaR,jet_pt,ee_js);

                        OptUn_eec->Fill(jet_pt,deltaR,ee_js_3D);
                        OptUn_eec->Fill(jet_pt,deltaR,ee_js_3D);

                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jss_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jss_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jss_3D);

                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jjs_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jjs_3D);
                        OptUn_e3c->Fill(jet_pt,deltaR,eee_jjs_3D);


                      
                }
            
           for(int k=s+1; k<int(fConstituents.size()) ; k++)
            { 
              // if(fCout)cout<<"in loop k"<<endl;

                float eee_jsk =((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[k].pt())/(pow(jet_pt,3)));

                float eee_jsk_3D =((fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[k].pt())/(pow(jet_pt,3)));
              
                float dR_js = delR(fConstituents[j],fConstituents[s]);
                float dR_sk = delR(fConstituents[s],fConstituents[k]);
                float dR_kj = delR(fConstituents[j],fConstituents[k]);

                float RL = -1;

                if(dR_js > dR_sk && dR_js > dR_kj){RL = dR_js;}
                else if(dR_sk > dR_js && dR_sk > dR_kj){RL = dR_sk;}
                else{RL = dR_kj;}

                 if(fpaircut == 1)
                {
                    float s_eta = fConstituents[s].eta();
                    float k_eta = fConstituents[k].eta();
                    float j_eta = fConstituents[j].eta();
                    float del_js_eta = abs(j_eta-s_eta);
                    float del_sk_eta = abs(s_eta-k_eta);
                    float del_kj_eta = abs(k_eta-j_eta);
                    if (del_sk_eta < 0.008 || del_js_eta < 0.008 || del_kj_eta < 0.008) continue;
                    else
                    {
                        E3C_pt_hist->Fill(RL,jet_pt,eee_jsk);

                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);

                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                    }
                }
                else
                {
                        E3C_pt_hist->Fill(RL,jet_pt,eee_jsk);

                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);

                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);
                        OptUn_e3c->Fill(jet_pt,RL,eee_jsk_3D);

                }
             }
        }
    }

 if(fCout)cout<<"Computed ENC"<<endl;
    
    //catch (fastjet::Error)
    // {
    //    AliError(" [w] FJ Exception caught.");
    //    // return -1;
    // } //end error message
    return;
}

//_________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::RunChanged(Int_t newrun)
{
  if(fStoreTrig)
  {
    auto downscalehandler = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    if(downscalehandler->GetCurrentRun() != newrun)
    {
      downscalehandler->SetRun(newrun);
    }
  }
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::RetrieveEventObjects()
{
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::Terminate(Option_t *)
{
  if(fCout) {std::cout << " Info::anrai: ===== In Terminate ===== " << std::endl;}
  // Called once at the end of the analysis.
  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}

