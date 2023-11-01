#include "TChain.h" 
#include "TTree.h"
#include "TList.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "TFile.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODTrackSelection.h"
#include "AliVAODHeader.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"

#include "AliKFParticleBase.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <string>

#include "TDatabasePDG.h"
#include "TPDGCode.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTask_Ld_CreateTrees_PairsOnly.h"

using namespace std;
ClassImp(AliAnalysisTask_Ld_CreateTrees_PairsOnly) 






AliAnalysisTask_Ld_CreateTrees_PairsOnly::AliAnalysisTask_Ld_CreateTrees_PairsOnly() : AliAnalysisTaskSE(),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(0),
  fUseOpenCuts(0),
  fIsMC(0),
  fSavePairsOnly(0),
  fSaveTree_Lambda(0),
  fLambda_px(0),
  fLambda_py(0),
  fLambda_pz(0),
  fLambda_px_Generated(0),
  fLambda_py_Generated(0),
  fLambda_pz_Generated(0),
  fLambda_Eta(0),
  fLambda_Phi(0),
  fLambda_TransverseRadius(0),
  fLambda_CosinePointingAngle(0),
  fLambda_DCAv0ToPrimaryVertex(0),
  fLambda_DCAv0Daughters(0),
  fLambda_Alpha(0),
  fLambda_qT(0),
  fLambda_DecayLength(0),
  fLambda_OpenAngle(0),
  fLambda_PDG_Daughter1(0),
  fLambda_PDG_Daughter2(0),
  fLambda_PDG_v01(0),
  fLambda_PDG_v02(0),
  fLambda_PDG_Mother1(0),
  fLambda_PDG_Mother2(0),
  fLambda_SameV0(0),
  fLambda_Event_Centrality(0),
  fLambda_Event_PrimaryVertexZ(0),
  fLambda_Event_BField(0),
  fLambda_Event_Multiplicity(0),
  fLambda_Event_Identifier(0),
  fLambda_Event_ParticleNumber(0),
  fLambda_Event_nParticles(0),
  fLambda_Daughter_Proton_px(0),
  fLambda_Daughter_Proton_py(0),
  fLambda_Daughter_Proton_pz(0),
  fLambda_Daughter_Proton_px_Generated(0),
  fLambda_Daughter_Proton_py_Generated(0),
  fLambda_Daughter_Proton_pz_Generated(0),
  fLambda_Daughter_Proton_px_DecayVertex(0),
  fLambda_Daughter_Proton_py_DecayVertex(0),
  fLambda_Daughter_Proton_pz_DecayVertex(0),
  fLambda_Daughter_Proton_pTPC(0),
  fLambda_Daughter_Proton_Eta(0),
  fLambda_Daughter_Proton_Phi(0),
  fLambda_Daughter_Proton_TPC_Chi2(0),
  fLambda_Daughter_Proton_TPC_dEdx(0),
  fLambda_Daughter_Proton_TPC_dEdx_nSigma(0),
  fLambda_Daughter_Proton_TOF_Mass2(0),
  fLambda_Daughter_Proton_TOF_Mass2_nSigma(0),
  fLambda_Daughter_Proton_ITS_dEdx(0),
  fLambda_Daughter_Proton_ITS_dEdx_nSigma(0),
  fLambda_Daughter_Proton_DCAxy(0),
  fLambda_Daughter_Proton_DCAz(0),
  fLambda_Daughter_Proton_TPC_nCrossedRows(0),
  fLambda_Daughter_Proton_TPC_nSharedCluster(0),
  fLambda_Daughter_Proton_TPC_nFindableCluster(0),
  fLambda_Daughter_Proton_TPC_nCluster(0),
  fLambda_Daughter_Proton_ITS_nCluster(0),
  fLambda_Daughter_Proton_ID(0),
  fLambda_Daughter_AntiPion_px(0),
  fLambda_Daughter_AntiPion_py(0),
  fLambda_Daughter_AntiPion_pz(0),
  fLambda_Daughter_AntiPion_px_Generated(0),
  fLambda_Daughter_AntiPion_py_Generated(0),
  fLambda_Daughter_AntiPion_pz_Generated(0),
  fLambda_Daughter_AntiPion_px_DecayVertex(0),
  fLambda_Daughter_AntiPion_py_DecayVertex(0),
  fLambda_Daughter_AntiPion_pz_DecayVertex(0),
  fLambda_Daughter_AntiPion_pTPC(0),
  fLambda_Daughter_AntiPion_Eta(0),
  fLambda_Daughter_AntiPion_Phi(0),
  fLambda_Daughter_AntiPion_TPC_Chi2(0),
  fLambda_Daughter_AntiPion_TPC_dEdx(0),
  fLambda_Daughter_AntiPion_TPC_dEdx_nSigma(0),
  fLambda_Daughter_AntiPion_TOF_Mass2(0),
  fLambda_Daughter_AntiPion_TOF_Mass2_nSigma(0),
  fLambda_Daughter_AntiPion_ITS_dEdx(0),
  fLambda_Daughter_AntiPion_ITS_dEdx_nSigma(0),
  fLambda_Daughter_AntiPion_DCAxy(0),
  fLambda_Daughter_AntiPion_DCAz(0),
  fLambda_Daughter_AntiPion_TPC_nCrossedRows(0),
  fLambda_Daughter_AntiPion_TPC_nSharedCluster(0),
  fLambda_Daughter_AntiPion_TPC_nFindableCluster(0),
  fLambda_Daughter_AntiPion_TPC_nCluster(0),
  fLambda_Daughter_AntiPion_ITS_nCluster(0),
  fLambda_Daughter_AntiPion_ID(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(0),
  fDeuteron_py(0),
  fDeuteron_pz(0),
  fDeuteron_px_Generated(0),
  fDeuteron_py_Generated(0),
  fDeuteron_pz_Generated(0),
  fDeuteron_pTPC(0),
  fDeuteron_Eta(0),
  fDeuteron_Phi(0),
  fDeuteron_TPC_Chi2(0),
  fDeuteron_TPC_dEdx(0),
  fDeuteron_TPC_dEdx_nSigma(0),
  fDeuteron_TOF_Mass2(0),
  fDeuteron_TOF_Mass2_nSigma(0),
  fDeuteron_ITS_dEdx(0),
  fDeuteron_ITS_dEdx_nSigma(0),
  fDeuteron_DCAxy(0),
  fDeuteron_DCAz(0),
  fDeuteron_Event_Centrality(0),
  fDeuteron_Event_PrimaryVertexZ(0),
  fDeuteron_Event_BField(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nSharedCluster(0),
  fDeuteron_TPC_nFindableCluster(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_PDG(0),
  fDeuteron_MotherPDG(0),
  fDeuteron_FilterBit(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Multiplicity(0),
  fDeuteron_Event_Identifier(0),
  fDeuteron_Event_ParticleNumber(0),
  fDeuteron_Event_nParticles(0),
  fSaveTree_AntiLambda(0),
  fAntiLambda_px(0),
  fAntiLambda_py(0),
  fAntiLambda_pz(0),
  fAntiLambda_px_Generated(0),
  fAntiLambda_py_Generated(0),
  fAntiLambda_pz_Generated(0),
  fAntiLambda_Eta(0),
  fAntiLambda_Phi(0),
  fAntiLambda_TransverseRadius(0),
  fAntiLambda_CosinePointingAngle(0),
  fAntiLambda_DCAv0ToPrimaryVertex(0),
  fAntiLambda_DCAv0Daughters(0),
  fAntiLambda_Alpha(0),
  fAntiLambda_qT(0),
  fAntiLambda_DecayLength(0),
  fAntiLambda_OpenAngle(0),
  fAntiLambda_PDG_Daughter1(0),
  fAntiLambda_PDG_Daughter2(0),
  fAntiLambda_PDG_v01(0),
  fAntiLambda_PDG_v02(0),
  fAntiLambda_PDG_Mother1(0),
  fAntiLambda_PDG_Mother2(0),
  fAntiLambda_SameV0(0),
  fAntiLambda_Event_Centrality(0),
  fAntiLambda_Event_PrimaryVertexZ(0),
  fAntiLambda_Event_BField(0),
  fAntiLambda_Event_Multiplicity(0),
  fAntiLambda_Event_Identifier(0),
  fAntiLambda_Event_ParticleNumber(0),
  fAntiLambda_Event_nParticles(0),
  fAntiLambda_Daughter_AntiProton_px(0),
  fAntiLambda_Daughter_AntiProton_py(0),
  fAntiLambda_Daughter_AntiProton_pz(0),
  fAntiLambda_Daughter_AntiProton_px_DecayVertex(0),
  fAntiLambda_Daughter_AntiProton_py_DecayVertex(0),
  fAntiLambda_Daughter_AntiProton_pz_DecayVertex(0),
  fAntiLambda_Daughter_AntiProton_pTPC(0),
  fAntiLambda_Daughter_AntiProton_Eta(0),
  fAntiLambda_Daughter_AntiProton_Phi(0),
  fAntiLambda_Daughter_AntiProton_TPC_Chi2(0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx(0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma(0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2(0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma(0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx(0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma(0),
  fAntiLambda_Daughter_AntiProton_DCAxy(0),
  fAntiLambda_Daughter_AntiProton_DCAz(0),
  fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_AntiProton_TPC_nSharedCluster(0),
  fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_AntiProton_TPC_nCluster(0),
  fAntiLambda_Daughter_AntiProton_ITS_nCluster(0),
  fAntiLambda_Daughter_AntiProton_ID(0),
  fAntiLambda_Daughter_Pion_px(0),
  fAntiLambda_Daughter_Pion_py(0),
  fAntiLambda_Daughter_Pion_pz(0),
  fAntiLambda_Daughter_Pion_px_DecayVertex(0),
  fAntiLambda_Daughter_Pion_py_DecayVertex(0),
  fAntiLambda_Daughter_Pion_pz_DecayVertex(0),
  fAntiLambda_Daughter_Pion_pTPC(0),
  fAntiLambda_Daughter_Pion_Eta(0),
  fAntiLambda_Daughter_Pion_Phi(0),
  fAntiLambda_Daughter_Pion_TPC_Chi2(0),
  fAntiLambda_Daughter_Pion_TPC_dEdx(0),
  fAntiLambda_Daughter_Pion_TPC_dEdx_nSigma(0),
  fAntiLambda_Daughter_Pion_TOF_Mass2(0),
  fAntiLambda_Daughter_Pion_TOF_Mass2_nSigma(0),
  fAntiLambda_Daughter_Pion_ITS_dEdx(0),
  fAntiLambda_Daughter_Pion_ITS_dEdx_nSigma(0),
  fAntiLambda_Daughter_Pion_DCAxy(0),
  fAntiLambda_Daughter_Pion_DCAz(0),
  fAntiLambda_Daughter_Pion_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_Pion_TPC_nSharedCluster(0),
  fAntiLambda_Daughter_Pion_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_Pion_TPC_nCluster(0),
  fAntiLambda_Daughter_Pion_ITS_nCluster(0),
  fAntiLambda_Daughter_Pion_ID(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(0),
  fAntiDeuteron_py(0),
  fAntiDeuteron_pz(0),
  fAntiDeuteron_px_Generated(0),
  fAntiDeuteron_py_Generated(0),
  fAntiDeuteron_pz_Generated(0),
  fAntiDeuteron_pTPC(0),
  fAntiDeuteron_Eta(0),
  fAntiDeuteron_Phi(0),
  fAntiDeuteron_TPC_Chi2(0),
  fAntiDeuteron_TPC_dEdx(0),
  fAntiDeuteron_TPC_dEdx_nSigma(0),
  fAntiDeuteron_TOF_Mass2(0),
  fAntiDeuteron_TOF_Mass2_nSigma(0),
  fAntiDeuteron_ITS_dEdx(0),
  fAntiDeuteron_ITS_dEdx_nSigma(0),
  fAntiDeuteron_DCAxy(0),
  fAntiDeuteron_DCAz(0),
  fAntiDeuteron_Event_Centrality(0),
  fAntiDeuteron_Event_PrimaryVertexZ(0),
  fAntiDeuteron_Event_BField(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nSharedCluster(0),
  fAntiDeuteron_TPC_nFindableCluster(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_PDG(0),
  fAntiDeuteron_MotherPDG(0),
  fAntiDeuteron_FilterBit(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Multiplicity(0),
  fAntiDeuteron_Event_Identifier(0),
  fAntiDeuteron_Event_ParticleNumber(0),
  fAntiDeuteron_Event_nParticles(0),
  fHistoList(0),
  h_Proton_TOF_m2_NoTOFcut(0),
  h_Deuteron_TOF_m2_NoTOFcut(0),
  h_AntiProton_TOF_m2_NoTOFcut(0),
  h_AntiDeuteron_TOF_m2_NoTOFcut(0),
  h_Proton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut(0),
  h_Pion_TOF_m2_NoTOFcut(0),
  h_AntiPion_TOF_m2_NoTOFcut(0),
  h_Pion_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiPion_ITS_dEdx_NoTOFcutNoITScut(0)
{


}



AliAnalysisTask_Ld_CreateTrees_PairsOnly::AliAnalysisTask_Ld_CreateTrees_PairsOnly(const char *name,int CollisionSystem, bool UseOpenCuts, bool isMC, bool SavePairsOnly) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(CollisionSystem),
  fUseOpenCuts(UseOpenCuts),
  fIsMC(isMC),
  fSavePairsOnly(SavePairsOnly),
  fSaveTree_Lambda(0),
  fLambda_px(0),
  fLambda_py(0),
  fLambda_pz(0),
  fLambda_px_Generated(0),
  fLambda_py_Generated(0),
  fLambda_pz_Generated(0),
  fLambda_Eta(0),
  fLambda_Phi(0),
  fLambda_TransverseRadius(0),
  fLambda_CosinePointingAngle(0),
  fLambda_DCAv0ToPrimaryVertex(0),
  fLambda_DCAv0Daughters(0),
  fLambda_Alpha(0),
  fLambda_qT(0),
  fLambda_DecayLength(0),
  fLambda_OpenAngle(0),
  fLambda_PDG_Daughter1(0),
  fLambda_PDG_Daughter2(0),
  fLambda_PDG_v01(0),
  fLambda_PDG_v02(0),
  fLambda_PDG_Mother1(0),
  fLambda_PDG_Mother2(0),
  fLambda_SameV0(0),
  fLambda_Event_Centrality(0),
  fLambda_Event_PrimaryVertexZ(0),
  fLambda_Event_BField(0),
  fLambda_Event_Multiplicity(0),
  fLambda_Event_Identifier(0),
  fLambda_Event_ParticleNumber(0),
  fLambda_Event_nParticles(0),
  fLambda_Daughter_Proton_px(0),
  fLambda_Daughter_Proton_py(0),
  fLambda_Daughter_Proton_pz(0),
  fLambda_Daughter_Proton_px_Generated(0),
  fLambda_Daughter_Proton_py_Generated(0),
  fLambda_Daughter_Proton_pz_Generated(0),
  fLambda_Daughter_Proton_px_DecayVertex(0),
  fLambda_Daughter_Proton_py_DecayVertex(0),
  fLambda_Daughter_Proton_pz_DecayVertex(0),
  fLambda_Daughter_Proton_pTPC(0),
  fLambda_Daughter_Proton_Eta(0),
  fLambda_Daughter_Proton_Phi(0),
  fLambda_Daughter_Proton_TPC_Chi2(0),
  fLambda_Daughter_Proton_TPC_dEdx(0),
  fLambda_Daughter_Proton_TPC_dEdx_nSigma(0),
  fLambda_Daughter_Proton_TOF_Mass2(0),
  fLambda_Daughter_Proton_TOF_Mass2_nSigma(0),
  fLambda_Daughter_Proton_ITS_dEdx(0),
  fLambda_Daughter_Proton_ITS_dEdx_nSigma(0),
  fLambda_Daughter_Proton_DCAxy(0),
  fLambda_Daughter_Proton_DCAz(0),
  fLambda_Daughter_Proton_TPC_nCrossedRows(0),
  fLambda_Daughter_Proton_TPC_nSharedCluster(0),
  fLambda_Daughter_Proton_TPC_nFindableCluster(0),
  fLambda_Daughter_Proton_TPC_nCluster(0),
  fLambda_Daughter_Proton_ITS_nCluster(0),
  fLambda_Daughter_Proton_ID(0),
  fLambda_Daughter_AntiPion_px(0),
  fLambda_Daughter_AntiPion_py(0),
  fLambda_Daughter_AntiPion_pz(0),
  fLambda_Daughter_AntiPion_px_Generated(0),
  fLambda_Daughter_AntiPion_py_Generated(0),
  fLambda_Daughter_AntiPion_pz_Generated(0),
  fLambda_Daughter_AntiPion_px_DecayVertex(0),
  fLambda_Daughter_AntiPion_py_DecayVertex(0),
  fLambda_Daughter_AntiPion_pz_DecayVertex(0),
  fLambda_Daughter_AntiPion_pTPC(0),
  fLambda_Daughter_AntiPion_Eta(0),
  fLambda_Daughter_AntiPion_Phi(0),
  fLambda_Daughter_AntiPion_TPC_Chi2(0),
  fLambda_Daughter_AntiPion_TPC_dEdx(0),
  fLambda_Daughter_AntiPion_TPC_dEdx_nSigma(0),
  fLambda_Daughter_AntiPion_TOF_Mass2(0),
  fLambda_Daughter_AntiPion_TOF_Mass2_nSigma(0),
  fLambda_Daughter_AntiPion_ITS_dEdx(0),
  fLambda_Daughter_AntiPion_ITS_dEdx_nSigma(0),
  fLambda_Daughter_AntiPion_DCAxy(0),
  fLambda_Daughter_AntiPion_DCAz(0),
  fLambda_Daughter_AntiPion_TPC_nCrossedRows(0),
  fLambda_Daughter_AntiPion_TPC_nSharedCluster(0),
  fLambda_Daughter_AntiPion_TPC_nFindableCluster(0),
  fLambda_Daughter_AntiPion_TPC_nCluster(0),
  fLambda_Daughter_AntiPion_ITS_nCluster(0),
  fLambda_Daughter_AntiPion_ID(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(0),
  fDeuteron_py(0),
  fDeuteron_pz(0),
  fDeuteron_px_Generated(0),
  fDeuteron_py_Generated(0),
  fDeuteron_pz_Generated(0),
  fDeuteron_pTPC(0),
  fDeuteron_Eta(0),
  fDeuteron_Phi(0),
  fDeuteron_TPC_Chi2(0),
  fDeuteron_TPC_dEdx(0),
  fDeuteron_TPC_dEdx_nSigma(0),
  fDeuteron_TOF_Mass2(0),
  fDeuteron_TOF_Mass2_nSigma(0),
  fDeuteron_ITS_dEdx(0),
  fDeuteron_ITS_dEdx_nSigma(0),
  fDeuteron_DCAxy(0),
  fDeuteron_DCAz(0),
  fDeuteron_Event_Centrality(0),
  fDeuteron_Event_PrimaryVertexZ(0),
  fDeuteron_Event_BField(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nSharedCluster(0),
  fDeuteron_TPC_nFindableCluster(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_PDG(0),
  fDeuteron_MotherPDG(0),
  fDeuteron_FilterBit(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Multiplicity(0),
  fDeuteron_Event_Identifier(0),
  fDeuteron_Event_ParticleNumber(0),
  fDeuteron_Event_nParticles(0),
  fSaveTree_AntiLambda(0),
  fAntiLambda_px(0),
  fAntiLambda_py(0),
  fAntiLambda_pz(0),
  fAntiLambda_px_Generated(0),
  fAntiLambda_py_Generated(0),
  fAntiLambda_pz_Generated(0),
  fAntiLambda_Eta(0),
  fAntiLambda_Phi(0),
  fAntiLambda_TransverseRadius(0),
  fAntiLambda_CosinePointingAngle(0),
  fAntiLambda_DCAv0ToPrimaryVertex(0),
  fAntiLambda_DCAv0Daughters(0),
  fAntiLambda_Alpha(0),
  fAntiLambda_qT(0),
  fAntiLambda_DecayLength(0),
  fAntiLambda_OpenAngle(0),
  fAntiLambda_PDG_Daughter1(0),
  fAntiLambda_PDG_Daughter2(0),
  fAntiLambda_PDG_v01(0),
  fAntiLambda_PDG_v02(0),
  fAntiLambda_PDG_Mother1(0),
  fAntiLambda_PDG_Mother2(0),
  fAntiLambda_SameV0(0),
  fAntiLambda_Event_Centrality(0),
  fAntiLambda_Event_PrimaryVertexZ(0),
  fAntiLambda_Event_BField(0),
  fAntiLambda_Event_Multiplicity(0),
  fAntiLambda_Event_Identifier(0),
  fAntiLambda_Event_ParticleNumber(0),
  fAntiLambda_Event_nParticles(0),
  fAntiLambda_Daughter_AntiProton_px(0),
  fAntiLambda_Daughter_AntiProton_py(0),
  fAntiLambda_Daughter_AntiProton_pz(0),
  fAntiLambda_Daughter_AntiProton_px_Generated(0),
  fAntiLambda_Daughter_AntiProton_py_Generated(0),
  fAntiLambda_Daughter_AntiProton_pz_Generated(0),
  fAntiLambda_Daughter_AntiProton_px_DecayVertex(0),
  fAntiLambda_Daughter_AntiProton_py_DecayVertex(0),
  fAntiLambda_Daughter_AntiProton_pz_DecayVertex(0),
  fAntiLambda_Daughter_AntiProton_pTPC(0),
  fAntiLambda_Daughter_AntiProton_Eta(0),
  fAntiLambda_Daughter_AntiProton_Phi(0),
  fAntiLambda_Daughter_AntiProton_TPC_Chi2(0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx(0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma(0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2(0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma(0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx(0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma(0),
  fAntiLambda_Daughter_AntiProton_DCAxy(0),
  fAntiLambda_Daughter_AntiProton_DCAz(0),
  fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_AntiProton_TPC_nSharedCluster(0),
  fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_AntiProton_TPC_nCluster(0),
  fAntiLambda_Daughter_AntiProton_ITS_nCluster(0),
  fAntiLambda_Daughter_AntiProton_ID(0),
  fAntiLambda_Daughter_Pion_px(0),
  fAntiLambda_Daughter_Pion_py(0),
  fAntiLambda_Daughter_Pion_pz(0),
  fAntiLambda_Daughter_Pion_px_Generated(0),
  fAntiLambda_Daughter_Pion_py_Generated(0),
  fAntiLambda_Daughter_Pion_pz_Generated(0),
  fAntiLambda_Daughter_Pion_px_DecayVertex(0),
  fAntiLambda_Daughter_Pion_py_DecayVertex(0),
  fAntiLambda_Daughter_Pion_pz_DecayVertex(0),
  fAntiLambda_Daughter_Pion_pTPC(0),
  fAntiLambda_Daughter_Pion_Eta(0),
  fAntiLambda_Daughter_Pion_Phi(0),
  fAntiLambda_Daughter_Pion_TPC_Chi2(0),
  fAntiLambda_Daughter_Pion_TPC_dEdx(0),
  fAntiLambda_Daughter_Pion_TPC_dEdx_nSigma(0),
  fAntiLambda_Daughter_Pion_TOF_Mass2(0),
  fAntiLambda_Daughter_Pion_TOF_Mass2_nSigma(0),
  fAntiLambda_Daughter_Pion_ITS_dEdx(0),
  fAntiLambda_Daughter_Pion_ITS_dEdx_nSigma(0),
  fAntiLambda_Daughter_Pion_DCAxy(0),
  fAntiLambda_Daughter_Pion_DCAz(0),
  fAntiLambda_Daughter_Pion_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_Pion_TPC_nSharedCluster(0),
  fAntiLambda_Daughter_Pion_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_Pion_TPC_nCluster(0),
  fAntiLambda_Daughter_Pion_ITS_nCluster(0),
  fAntiLambda_Daughter_Pion_ID(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(0),
  fAntiDeuteron_py(0),
  fAntiDeuteron_pz(0),
  fAntiDeuteron_px_Generated(0),
  fAntiDeuteron_py_Generated(0),
  fAntiDeuteron_pz_Generated(0),
  fAntiDeuteron_pTPC(0),
  fAntiDeuteron_Eta(0),
  fAntiDeuteron_Phi(0),
  fAntiDeuteron_TPC_Chi2(0),
  fAntiDeuteron_TPC_dEdx(0),
  fAntiDeuteron_TPC_dEdx_nSigma(0),
  fAntiDeuteron_TOF_Mass2(0),
  fAntiDeuteron_TOF_Mass2_nSigma(0),
  fAntiDeuteron_ITS_dEdx(0),
  fAntiDeuteron_ITS_dEdx_nSigma(0),
  fAntiDeuteron_DCAxy(0),
  fAntiDeuteron_DCAz(0),
  fAntiDeuteron_Event_Centrality(0),
  fAntiDeuteron_Event_PrimaryVertexZ(0),
  fAntiDeuteron_Event_BField(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nSharedCluster(0),
  fAntiDeuteron_TPC_nFindableCluster(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_PDG(0),
  fAntiDeuteron_MotherPDG(0),
  fAntiDeuteron_FilterBit(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Multiplicity(0),
  fAntiDeuteron_Event_Identifier(0),
  fAntiDeuteron_Event_ParticleNumber(0),
  fAntiDeuteron_Event_nParticles(0),
  fHistoList(0),
  h_Proton_TOF_m2_NoTOFcut(0),
  h_Deuteron_TOF_m2_NoTOFcut(0),
  h_AntiProton_TOF_m2_NoTOFcut(0),
  h_AntiDeuteron_TOF_m2_NoTOFcut(0),
  h_Proton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut(0),
  h_Pion_TOF_m2_NoTOFcut(0),
  h_AntiPion_TOF_m2_NoTOFcut(0),
  h_Pion_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiPion_ITS_dEdx_NoTOFcutNoITScut(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());
  DefineOutput(5,TList::Class());

}

  
AliAnalysisTask_Ld_CreateTrees_PairsOnly::~AliAnalysisTask_Ld_CreateTrees_PairsOnly()
{

  if(fSaveTree_Lambda)
    {
      delete fSaveTree_Lambda;
    }

  if(fSaveTree_Deuteron)
    {
      delete fSaveTree_Deuteron;
    }

  if(fSaveTree_AntiLambda)
    {
      delete fSaveTree_AntiLambda;
    }

  if(fSaveTree_AntiDeuteron)
    {
      delete fSaveTree_AntiDeuteron;
    }

  if(fHistoList)
    {
      delete fHistoList;
    }

}







void AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserCreateOutputObjects()
{

  fHistoList = new TList();
  fHistoList->SetOwner();

  h_Proton_TOF_m2_NoTOFcut = new TH2F("h_Proton_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (protons)",240,0.0,6.0,500,0.0,10.0);
  h_Proton_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_Proton_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_Proton_TOF_m2_NoTOFcut);

  h_Deuteron_TOF_m2_NoTOFcut = new TH2F("h_Deuteron_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (deuterons)",240,0.0,6.0,500,0.0,10.0);
  h_Deuteron_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_Deuteron_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_Deuteron_TOF_m2_NoTOFcut);

  h_AntiProton_TOF_m2_NoTOFcut = new TH2F("h_AntiProton_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (antiprotons)",240,0.0,6.0,500,0.0,10.0);
  h_AntiProton_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_AntiProton_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_AntiProton_TOF_m2_NoTOFcut);

  h_AntiDeuteron_TOF_m2_NoTOFcut = new TH2F("h_AntiDeuteron_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (antideuterons)",240,0.0,6.0,500,0.0,10.0);
  h_AntiDeuteron_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_AntiDeuteron_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_AntiDeuteron_TOF_m2_NoTOFcut);


  h_Proton_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_Proton_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (protons)",400,0.0,10.0,500,0.0,500.0);
  h_Proton_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_Proton_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_Proton_ITS_dEdx_NoTOFcutNoITScut);

  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_Deuteron_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (deuterons)",400,0.0,10.0,500,0.0,500.0);
  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_Deuteron_ITS_dEdx_NoTOFcutNoITScut);

  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_AntiProton_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (antiprotons)",400,0.0,10.0,500,0.0,500.0);
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_AntiProton_ITS_dEdx_NoTOFcutNoITScut);

  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (antideuterons)",400,0.0,10.0,500,0.0,500.0);
  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut);

  h_Pion_TOF_m2_NoTOFcut = new TH2F("h_Pion_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (Pions)",240,0.0,6.0,600,0.0,0.3);
  h_Pion_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_Pion_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_Pion_TOF_m2_NoTOFcut);

  h_AntiPion_TOF_m2_NoTOFcut = new TH2F("h_AntiPion_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (AntiPions)",240,0.0,6.0,600,0.0,0.3);
  h_AntiPion_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_AntiPion_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_AntiPion_TOF_m2_NoTOFcut);

  h_Pion_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_Pion_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (Pions)",400,0.0,10.0,500,0.0,500.0);
  h_Pion_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_Pion_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_Pion_ITS_dEdx_NoTOFcutNoITScut);

  h_AntiPion_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_AntiPion_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (AntiPions)",400,0.0,10.0,500,0.0,500.0);
  h_AntiPion_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_AntiPion_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_AntiPion_ITS_dEdx_NoTOFcutNoITScut);


  fSaveTree_Lambda = new TTree("fSaveTree_Lambda","fSaveTree_Lambda");
  fSaveTree_Lambda->Branch("Lambda_px",&fLambda_px,"Lambda_px/F");
  fSaveTree_Lambda->Branch("Lambda_py",&fLambda_py,"Lambda_py/F");
  fSaveTree_Lambda->Branch("Lambda_pz",&fLambda_pz,"Lambda_pz/F");
  fSaveTree_Lambda->Branch("Lambda_px_Generated",&fLambda_px_Generated,"Lambda_px_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_py_Generated",&fLambda_py_Generated,"Lambda_py_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_pz_Generated",&fLambda_pz_Generated,"Lambda_pz_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Eta",&fLambda_Eta,"Lambda_Eta/F");
  fSaveTree_Lambda->Branch("Lambda_Phi",&fLambda_Phi,"Lambda_Phi/F");
  fSaveTree_Lambda->Branch("Lambda_TransverseRadius",&fLambda_TransverseRadius,"Lambda_TransverseRadius/F");
  fSaveTree_Lambda->Branch("Lambda_CosinePointingAngle",&fLambda_CosinePointingAngle,"Lambda_CosinePointingAngle/F");
  fSaveTree_Lambda->Branch("Lambda_DCAv0ToPrimaryVertex",&fLambda_DCAv0ToPrimaryVertex,"Lambda_DCAv0ToPrimaryVertex/F");
  fSaveTree_Lambda->Branch("Lambda_DCAv0Daughters",&fLambda_DCAv0Daughters,"Lambda_DCAv0Daughters/F");
  fSaveTree_Lambda->Branch("Lambda_Alpha",&fLambda_Alpha,"Lambda_Alpha/F");
  fSaveTree_Lambda->Branch("Lambda_qT",&fLambda_qT,"Lambda_qT/F");
  fSaveTree_Lambda->Branch("Lambda_DecayLength",&fLambda_DecayLength,"Lambda_DecayLength/F");
  fSaveTree_Lambda->Branch("Lambda_OpenAngle",&fLambda_OpenAngle,"Lambda_OpenAngle/F");
  fSaveTree_Lambda->Branch("Lambda_PDG_Daughter1",&fLambda_PDG_Daughter1,"Lambda_PDG_Daughter1/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_Daughter2",&fLambda_PDG_Daughter2,"Lambda_PDG_Daughter2/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_v01",&fLambda_PDG_v01,"Lambda_PDG_v01/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_v02",&fLambda_PDG_v02,"Lambda_PDG_v02/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_Mother1",&fLambda_PDG_Mother1,"Lambda_PDG_Mother1/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_Mother2",&fLambda_PDG_Mother2,"Lambda_PDG_Mother2/I");
  fSaveTree_Lambda->Branch("Lambda_SameV0",&fLambda_SameV0,"Lambda_SameV0/O");
  fSaveTree_Lambda->Branch("Lambda_Event_Centrality",&fLambda_Event_Centrality,"Lambda_Event_Centrality/F");
  fSaveTree_Lambda->Branch("Lambda_Event_PrimaryVertexZ",&fLambda_Event_PrimaryVertexZ,"Lambda_Event_PrimaryVertexZ/F");
  fSaveTree_Lambda->Branch("Lambda_Event_BField",&fLambda_Event_BField,"Lambda_Event_BField/F");
  fSaveTree_Lambda->Branch("Lambda_Event_Multiplicity",&fLambda_Event_Multiplicity,"Lambda_Event_Multiplicity/i");
  fSaveTree_Lambda->Branch("Lambda_Event_Identifier",&fLambda_Event_Identifier,"Lambda_Event_Identifier/l");
  fSaveTree_Lambda->Branch("Lambda_Event_ParticleNumber",&fLambda_Event_ParticleNumber,"Lambda_Event_ParticleNumber/s");
  fSaveTree_Lambda->Branch("Lambda_Event_nParticles",&fLambda_Event_nParticles,"Lambda_Event_nParticles/s");

  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_px",&fLambda_Daughter_Proton_px,"Lambda_Daughter_Proton_px/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_py",&fLambda_Daughter_Proton_py,"Lambda_Daughter_Proton_py/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pz",&fLambda_Daughter_Proton_pz,"Lambda_Daughter_Proton_pz/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_px_Generated",&fLambda_Daughter_Proton_px_Generated,"Lambda_Daughter_Proton_px_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_py_Generated",&fLambda_Daughter_Proton_py_Generated,"Lambda_Daughter_Proton_py_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pz_Generated",&fLambda_Daughter_Proton_pz_Generated,"Lambda_Daughter_Proton_pz_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_px_DecayVertex",&fLambda_Daughter_Proton_px_DecayVertex,"Lambda_Daughter_Proton_px_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_py_DecayVertex",&fLambda_Daughter_Proton_py_DecayVertex,"Lambda_Daughter_Proton_py_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pz_DecayVertex",&fLambda_Daughter_Proton_pz_DecayVertex,"Lambda_Daughter_Proton_pz_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pTPC",&fLambda_Daughter_Proton_pTPC,"Lambda_Daughter_Proton_pTPC/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_Eta",&fLambda_Daughter_Proton_Eta,"Lambda_Daughter_Proton_Eta/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_Phi",&fLambda_Daughter_Proton_Phi,"Lambda_Daughter_Proton_Phi/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_Chi2",&fLambda_Daughter_Proton_TPC_Chi2,"Lambda_Daughter_Proton_TPC_Chi2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_dEdx",&fLambda_Daughter_Proton_TPC_dEdx,"Lambda_Daughter_Proton_TPC_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_dEdx_nSigma",&fLambda_Daughter_Proton_TPC_dEdx_nSigma,"Lambda_Daughter_Proton_TPC_dEdx_nSigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TOF_Mass2",&fLambda_Daughter_Proton_TOF_Mass2,"Lambda_Daughter_Proton_TOF_Mass2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TOF_Mass2_nSigma",&fLambda_Daughter_Proton_TOF_Mass2_nSigma,"Lambda_Daughter_Proton_TOF_Mass2_nSigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_dEdx",&fLambda_Daughter_Proton_ITS_dEdx,"Lambda_Daughter_Proton_ITS_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_dEdx_nSigma",&fLambda_Daughter_Proton_ITS_dEdx_nSigma,"Lambda_Daughter_Proton_ITS_dEdx_nSigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_DCAxy",&fLambda_Daughter_Proton_DCAxy,"Lambda_Daughter_Proton_DCAxy/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_DCAz",&fLambda_Daughter_Proton_DCAz,"Lambda_Daughter_Proton_DCAz/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nCrossedRows",&fLambda_Daughter_Proton_TPC_nCrossedRows,"Lambda_Daughter_Proton_TPC_nCrossedRows/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nSharedCluster",&fLambda_Daughter_Proton_TPC_nSharedCluster,"Lambda_Daughter_Proton_TPC_nSharedCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nFindableCluster",&fLambda_Daughter_Proton_TPC_nFindableCluster,"Lambda_Daughter_Proton_TPC_nFindableCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nCluster",&fLambda_Daughter_Proton_TPC_nCluster,"Lambda_Daughter_Proton_TPC_nCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_nCluster",&fLambda_Daughter_Proton_ITS_nCluster,"Lambda_Daughter_Proton_ITS_nCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ID",&fLambda_Daughter_Proton_ID,"Lambda_Daughter_Proton_ID/i");

  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_px",&fLambda_Daughter_AntiPion_px,"Lambda_Daughter_AntiPion_px/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_py",&fLambda_Daughter_AntiPion_py,"Lambda_Daughter_AntiPion_py/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz",&fLambda_Daughter_AntiPion_pz,"Lambda_Daughter_AntiPion_pz/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_px_Generated",&fLambda_Daughter_AntiPion_px_Generated,"Lambda_Daughter_AntiPion_px_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_py_Generated",&fLambda_Daughter_AntiPion_py_Generated,"Lambda_Daughter_AntiPion_py_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz_Generated",&fLambda_Daughter_AntiPion_pz_Generated,"Lambda_Daughter_AntiPion_pz_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_px_DecayVertex",&fLambda_Daughter_AntiPion_px_DecayVertex,"Lambda_Daughter_AntiPion_px_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_py_DecayVertex",&fLambda_Daughter_AntiPion_py_DecayVertex,"Lambda_Daughter_AntiPion_py_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz_DecayVertex",&fLambda_Daughter_AntiPion_pz_DecayVertex,"Lambda_Daughter_AntiPion_pz_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pTPC",&fLambda_Daughter_AntiPion_pTPC,"Lambda_Daughter_AntiPion_pTPC/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_Eta",&fLambda_Daughter_AntiPion_Eta,"Lambda_Daughter_AntiPion_Eta/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_Phi",&fLambda_Daughter_AntiPion_Phi,"Lambda_Daughter_AntiPion_Phi/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_Chi2",&fLambda_Daughter_AntiPion_TPC_Chi2,"Lambda_Daughter_AntiPion_TPC_Chi2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_dEdx",&fLambda_Daughter_AntiPion_TPC_dEdx,"Lambda_Daughter_AntiPion_TPC_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_dEdx_nSigma",&fLambda_Daughter_AntiPion_TPC_dEdx_nSigma,"Lambda_Daughter_AntiPion_TPC_dEdx_nSigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TOF_Mass2",&fLambda_Daughter_AntiPion_TOF_Mass2,"Lambda_Daughter_AntiPion_TOF_Mass2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TOF_Mass2_nSigma",&fLambda_Daughter_AntiPion_TOF_Mass2_nSigma,"Lambda_Daughter_AntiPion_TOF_Mass2_nSigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_dEdx",&fLambda_Daughter_AntiPion_ITS_dEdx,"Lambda_Daughter_AntiPion_ITS_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_dEdx_nSigma",&fLambda_Daughter_AntiPion_ITS_dEdx_nSigma,"Lambda_Daughter_AntiPion_ITS_dEdx_nSigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_DCAxy",&fLambda_Daughter_AntiPion_DCAxy,"Lambda_Daughter_AntiPion_DCAxy/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_DCAz",&fLambda_Daughter_AntiPion_DCAz,"Lambda_Daughter_AntiPion_DCAz/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nCrossedRows",&fLambda_Daughter_AntiPion_TPC_nCrossedRows,"Lambda_Daughter_AntiPion_TPC_nCrossedRows/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nSharedCluster",&fLambda_Daughter_AntiPion_TPC_nSharedCluster,"Lambda_Daughter_AntiPion_TPC_nSharedCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nFindableCluster",&fLambda_Daughter_AntiPion_TPC_nFindableCluster,"Lambda_Daughter_AntiPion_TPC_nFindableCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nCluster",&fLambda_Daughter_AntiPion_TPC_nCluster,"Lambda_Daughter_AntiPion_TPC_nCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_nCluster",&fLambda_Daughter_AntiPion_ITS_nCluster,"Lambda_Daughter_AntiPion_ITS_nCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ID",&fLambda_Daughter_AntiPion_ID,"Lambda_Daughter_AntiPion_ID/i");









  fSaveTree_Deuteron = new TTree("fSaveTree_Deuteron","fSaveTree_Deuteron");
  fSaveTree_Deuteron->Branch("Deuteron_px",&fDeuteron_px,"Deuteron_px/F");
  fSaveTree_Deuteron->Branch("Deuteron_py",&fDeuteron_py,"Deuteron_py/F");
  fSaveTree_Deuteron->Branch("Deuteron_pz",&fDeuteron_pz,"Deuteron_pz/F");
  fSaveTree_Deuteron->Branch("Deuteron_px_Generated",&fDeuteron_px_Generated,"Deuteron_px_Generated/F");
  fSaveTree_Deuteron->Branch("Deuteron_py_Generated",&fDeuteron_py_Generated,"Deuteron_py_Generated/F");
  fSaveTree_Deuteron->Branch("Deuteron_pz_Generated",&fDeuteron_pz_Generated,"Deuteron_pz_Generated/F");
  fSaveTree_Deuteron->Branch("Deuteron_pTPC",&fDeuteron_pTPC,"Deuteron_pTPC/F");
  fSaveTree_Deuteron->Branch("Deuteron_Eta",&fDeuteron_Eta,"Deuteron_Eta/F");
  fSaveTree_Deuteron->Branch("Deuteron_Phi",&fDeuteron_Phi,"Deuteron_Phi/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_Chi2",&fDeuteron_TPC_Chi2,"Deuteron_TPC_Chi2/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx",&fDeuteron_TPC_dEdx,"Deuteron_TPC_dEdx/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx_nSigma",&fDeuteron_TPC_dEdx_nSigma,"Deuteron_TPC_dEdx_nSigma/F");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2",&fDeuteron_TOF_Mass2,"Deuteron_TOF_Mass2/F");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2_nSigma",&fDeuteron_TOF_Mass2_nSigma,"Deuteron_TOF_Mass2_nSigma/F");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx",&fDeuteron_ITS_dEdx,"Deuteron_ITS_dEdx/F");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx_nSigma",&fDeuteron_ITS_dEdx_nSigma,"Deuteron_ITS_dEdx_nSigma/F");
  fSaveTree_Deuteron->Branch("Deuteron_DCAxy",&fDeuteron_DCAxy,"Deuteron_DCAxy/F");
  fSaveTree_Deuteron->Branch("Deuteron_DCAz",&fDeuteron_DCAz,"Deuteron_DCAz/F");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Centrality",&fDeuteron_Event_Centrality,"Deuteron_Event_Centrality/F");
  fSaveTree_Deuteron->Branch("Deuteron_Event_PrimaryVertexZ",&fDeuteron_Event_PrimaryVertexZ,"Deuteron_Event_PrimaryVertexZ/F");
  fSaveTree_Deuteron->Branch("Deuteron_Event_BField",&fDeuteron_Event_BField,"Deuteron_Event_BField/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&fDeuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&fDeuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nFindableCluster",&fDeuteron_TPC_nFindableCluster,"Deuteron_TPC_nFindableCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCluster",&fDeuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_nCluster",&fDeuteron_ITS_nCluster,"Deuteron_ITS_nCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_PDG",&fDeuteron_PDG,"Deuteron_PDG/I");
  fSaveTree_Deuteron->Branch("Deuteron_MotherPDG",&fDeuteron_MotherPDG,"Deuteron_MotherPDG/I");
  fSaveTree_Deuteron->Branch("Deuteron_FilterBit",&fDeuteron_FilterBit,"Deuteron_FilterBit/O");
  fSaveTree_Deuteron->Branch("Deuteron_ID",&fDeuteron_ID,"Deuteron_ID/i");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Multiplicity",&fDeuteron_Event_Multiplicity,"Deuteron_Event_Multiplicity/i");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Identifier",&fDeuteron_Event_Identifier,"Deuteron_Event_Identifier/l");
  fSaveTree_Deuteron->Branch("Deuteron_Event_ParticleNumber",&fDeuteron_Event_ParticleNumber,"Deuteron_Event_ParticleNumber/s");
  fSaveTree_Deuteron->Branch("Deuteron_Event_nParticles",&fDeuteron_Event_nParticles,"Deuteron_Event_nParticles/s");


  fSaveTree_AntiLambda = new TTree("fSaveTree_AntiLambda","fSaveTree_AntiLambda");
  fSaveTree_AntiLambda->Branch("AntiLambda_px",&fAntiLambda_px,"AntiLambda_px/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_py",&fAntiLambda_py,"AntiLambda_py/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_pz",&fAntiLambda_pz,"AntiLambda_pz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_px_Generated",&fAntiLambda_px_Generated,"AntiLambda_px_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_py_Generated",&fAntiLambda_py_Generated,"AntiLambda_py_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_pz_Generated",&fAntiLambda_pz_Generated,"AntiLambda_pz_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Eta",&fAntiLambda_Eta,"AntiLambda_Eta/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Phi",&fAntiLambda_Phi,"AntiLambda_Phi/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_TransverseRadius",&fAntiLambda_TransverseRadius,"AntiLambda_TransverseRadius/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_CosinePointingAngle",&fAntiLambda_CosinePointingAngle,"AntiLambda_CosinePointingAngle/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_DCAv0ToPrimaryVertex",&fAntiLambda_DCAv0ToPrimaryVertex,"AntiLambda_DCAv0ToPrimaryVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_DCAv0Daughters",&fAntiLambda_DCAv0Daughters,"AntiLambda_DCAv0Daughters/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Alpha",&fAntiLambda_Alpha,"AntiLambda_Alpha/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_qT",&fAntiLambda_qT,"AntiLambda_qT/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_DecayLength",&fAntiLambda_DecayLength,"AntiLambda_DecayLength/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_OpenAngle",&fAntiLambda_OpenAngle,"AntiLambda_OpenAngle/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Daughter1",&fAntiLambda_PDG_Daughter1,"AntiLambda_PDG_Daughter1/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Daughter2",&fAntiLambda_PDG_Daughter2,"AntiLambda_PDG_Daughter2/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_v01",&fAntiLambda_PDG_v01,"AntiLambda_PDG_v01/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_v02",&fAntiLambda_PDG_v02,"AntiLambda_PDG_v02/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Mother1",&fAntiLambda_PDG_Mother1,"AntiLambda_PDG_Mother1/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Mother2",&fAntiLambda_PDG_Mother2,"AntiLambda_PDG_Mother2/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_SameV0",&fAntiLambda_SameV0,"AntiLambda_SameV0/O");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_Centrality",&fAntiLambda_Event_Centrality,"AntiLambda_Event_Centrality/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_PrimaryVertexZ",&fAntiLambda_Event_PrimaryVertexZ,"AntiLambda_Event_PrimaryVertexZ/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_BField",&fAntiLambda_Event_BField,"AntiLambda_Event_BField/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_Multiplicity",&fAntiLambda_Event_Multiplicity,"AntiLambda_Event_Multiplicity/i");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_Identifier",&fAntiLambda_Event_Identifier,"AntiLambda_Event_Identifier/l");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_ParticleNumber",&fAntiLambda_Event_ParticleNumber,"AntiLambda_Event_ParticleNumber/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_nParticles",&fAntiLambda_Event_nParticles,"AntiLambda_Event_nParticles/s");


  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px",&fAntiLambda_Daughter_AntiProton_px,"AntiLambda_Daughter_AntiProton_px/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py",&fAntiLambda_Daughter_AntiProton_py,"AntiLambda_Daughter_AntiProton_py/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz",&fAntiLambda_Daughter_AntiProton_pz,"AntiLambda_Daughter_AntiProton_pz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px_Generated",&fAntiLambda_Daughter_AntiProton_px_Generated,"AntiLambda_Daughter_AntiProton_px_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py_Generated",&fAntiLambda_Daughter_AntiProton_py_Generated,"AntiLambda_Daughter_AntiProton_py_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz_Generated",&fAntiLambda_Daughter_AntiProton_pz_Generated,"AntiLambda_Daughter_AntiProton_pz_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px_DecayVertex",&fAntiLambda_Daughter_AntiProton_px_DecayVertex,"AntiLambda_Daughter_AntiProton_px_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py_DecayVertex",&fAntiLambda_Daughter_AntiProton_py_DecayVertex,"AntiLambda_Daughter_AntiProton_py_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz_DecayVertex",&fAntiLambda_Daughter_AntiProton_pz_DecayVertex,"AntiLambda_Daughter_AntiProton_pz_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pTPC",&fAntiLambda_Daughter_AntiProton_pTPC,"AntiLambda_Daughter_AntiProton_pTPC/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_Eta",&fAntiLambda_Daughter_AntiProton_Eta,"AntiLambda_Daughter_AntiProton_Eta/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_Phi",&fAntiLambda_Daughter_AntiProton_Phi,"AntiLambda_Daughter_AntiProton_Phi/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_Chi2",&fAntiLambda_Daughter_AntiProton_TPC_Chi2,"AntiLambda_Daughter_AntiProton_TPC_Chi2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_dEdx",&fAntiLambda_Daughter_AntiProton_TPC_dEdx,"AntiLambda_Daughter_AntiProton_TPC_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma",&fAntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma,"AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TOF_Mass2",&fAntiLambda_Daughter_AntiProton_TOF_Mass2,"AntiLambda_Daughter_AntiProton_TOF_Mass2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma",&fAntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma,"AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_dEdx",&fAntiLambda_Daughter_AntiProton_ITS_dEdx,"AntiLambda_Daughter_AntiProton_ITS_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma",&fAntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma,"AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_DCAxy",&fAntiLambda_Daughter_AntiProton_DCAxy,"AntiLambda_Daughter_AntiProton_DCAxy/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_DCAz",&fAntiLambda_Daughter_AntiProton_DCAz,"AntiLambda_Daughter_AntiProton_DCAz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nCrossedRows",&fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows,"AntiLambda_Daughter_AntiProton_TPC_nCrossedRows/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nSharedCluster",&fAntiLambda_Daughter_AntiProton_TPC_nSharedCluster,"AntiLambda_Daughter_AntiProton_TPC_nSharedCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nFindableCluster",&fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster,"AntiLambda_Daughter_AntiProton_TPC_nFindableCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nCluster",&fAntiLambda_Daughter_AntiProton_TPC_nCluster,"AntiLambda_Daughter_AntiProton_TPC_nCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_nCluster",&fAntiLambda_Daughter_AntiProton_ITS_nCluster,"AntiLambda_Daughter_AntiProton_ITS_nCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ID",&fAntiLambda_Daughter_AntiProton_ID,"AntiLambda_Daughter_AntiProton_ID/i");

  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px",&fAntiLambda_Daughter_Pion_px,"AntiLambda_Daughter_Pion_px/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py",&fAntiLambda_Daughter_Pion_py,"AntiLambda_Daughter_Pion_py/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz",&fAntiLambda_Daughter_Pion_pz,"AntiLambda_Daughter_Pion_pz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px_Generated",&fAntiLambda_Daughter_Pion_px_Generated,"AntiLambda_Daughter_Pion_px_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py_Generated",&fAntiLambda_Daughter_Pion_py_Generated,"AntiLambda_Daughter_Pion_py_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz_Generated",&fAntiLambda_Daughter_Pion_pz_Generated,"AntiLambda_Daughter_Pion_pz_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px_DecayVertex",&fAntiLambda_Daughter_Pion_px_DecayVertex,"AntiLambda_Daughter_Pion_px_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py_DecayVertex",&fAntiLambda_Daughter_Pion_py_DecayVertex,"AntiLambda_Daughter_Pion_py_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz_DecayVertex",&fAntiLambda_Daughter_Pion_pz_DecayVertex,"AntiLambda_Daughter_Pion_pz_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pTPC",&fAntiLambda_Daughter_Pion_pTPC,"AntiLambda_Daughter_Pion_pTPC/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_Eta",&fAntiLambda_Daughter_Pion_Eta,"AntiLambda_Daughter_Pion_Eta/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_Phi",&fAntiLambda_Daughter_Pion_Phi,"AntiLambda_Daughter_Pion_Phi/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_Chi2",&fAntiLambda_Daughter_Pion_TPC_Chi2,"AntiLambda_Daughter_Pion_TPC_Chi2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_dEdx",&fAntiLambda_Daughter_Pion_TPC_dEdx,"AntiLambda_Daughter_Pion_TPC_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_dEdx_nSigma",&fAntiLambda_Daughter_Pion_TPC_dEdx_nSigma,"AntiLambda_Daughter_Pion_TPC_dEdx_nSigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TOF_Mass2",&fAntiLambda_Daughter_Pion_TOF_Mass2,"AntiLambda_Daughter_Pion_TOF_Mass2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TOF_Mass2_nSigma",&fAntiLambda_Daughter_Pion_TOF_Mass2_nSigma,"AntiLambda_Daughter_Pion_TOF_Mass2_nSigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_dEdx",&fAntiLambda_Daughter_Pion_ITS_dEdx,"AntiLambda_Daughter_Pion_ITS_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_dEdx_nSigma",&fAntiLambda_Daughter_Pion_ITS_dEdx_nSigma,"AntiLambda_Daughter_Pion_ITS_dEdx_nSigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_DCAxy",&fAntiLambda_Daughter_Pion_DCAxy,"AntiLambda_Daughter_Pion_DCAxy/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_DCAz",&fAntiLambda_Daughter_Pion_DCAz,"AntiLambda_Daughter_Pion_DCAz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nCrossedRows",&fAntiLambda_Daughter_Pion_TPC_nCrossedRows,"AntiLambda_Daughter_Pion_TPC_nCrossedRows/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nSharedCluster",&fAntiLambda_Daughter_Pion_TPC_nSharedCluster,"AntiLambda_Daughter_Pion_TPC_nSharedCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nFindableCluster",&fAntiLambda_Daughter_Pion_TPC_nFindableCluster,"AntiLambda_Daughter_Pion_TPC_nFindableCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nCluster",&fAntiLambda_Daughter_Pion_TPC_nCluster,"AntiLambda_Daughter_Pion_TPC_nCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_nCluster",&fAntiLambda_Daughter_Pion_ITS_nCluster,"AntiLambda_Daughter_Pion_ITS_nCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ID",&fAntiLambda_Daughter_Pion_ID,"AntiLambda_Daughter_Pion_ID/i");



  fSaveTree_AntiDeuteron = new TTree("fSaveTree_AntiDeuteron","fSaveTree_AntiDeuteron");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_px",&fAntiDeuteron_px,"AntiDeuteron_px/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_py",&fAntiDeuteron_py,"AntiDeuteron_py/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pz",&fAntiDeuteron_pz,"AntiDeuteron_pz/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_px_Generated",&fAntiDeuteron_px_Generated,"AntiDeuteron_px_Generated/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_py_Generated",&fAntiDeuteron_py_Generated,"AntiDeuteron_py_Generated/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pz_Generated",&fAntiDeuteron_pz_Generated,"AntiDeuteron_pz_Generated/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&fAntiDeuteron_pTPC,"AntiDeuteron_pTPC/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&fAntiDeuteron_Eta,"AntiDeuteron_Eta/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&fAntiDeuteron_Phi,"AntiDeuteron_Phi/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2",&fAntiDeuteron_TPC_Chi2,"AntiDeuteron_TPC_Chi2/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&fAntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_nSigma",&fAntiDeuteron_TPC_dEdx_nSigma,"AntiDeuteron_TPC_dEdx_nSigma/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&fAntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_nSigma",&fAntiDeuteron_TOF_Mass2_nSigma,"AntiDeuteron_TOF_Mass2_nSigma/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&fAntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_nSigma",&fAntiDeuteron_ITS_dEdx_nSigma,"AntiDeuteron_ITS_dEdx_nSigma/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&fAntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&fAntiDeuteron_DCAz,"AntiDeuteron_DCAz/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Centrality",&fAntiDeuteron_Event_Centrality,"AntiDeuteron_Event_Centrality/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_PrimaryVertexZ",&fAntiDeuteron_Event_PrimaryVertexZ,"AntiDeuteron_Event_PrimaryVertexZ/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_BField",&fAntiDeuteron_Event_BField,"AntiDeuteron_Event_BField/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&fAntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&fAntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nFindableCluster",&fAntiDeuteron_TPC_nFindableCluster,"AntiDeuteron_TPC_nFindableCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&fAntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&fAntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_PDG",&fAntiDeuteron_PDG,"AntiDeuteron_PDG/I");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_MotherPDG",&fAntiDeuteron_MotherPDG,"AntiDeuteron_MotherPDG/I");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_FilterBit",&fAntiDeuteron_FilterBit,"AntiDeuteron_FilterBit/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ID",&fAntiDeuteron_ID,"AntiDeuteron_ID/i");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Multiplicity",&fAntiDeuteron_Event_Multiplicity,"AntiDeuteron_Event_Multiplicity/i");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&fAntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/l");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_ParticleNumber",&fAntiDeuteron_Event_ParticleNumber,"AntiDeuteron_Event_ParticleNumber/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_nParticles",&fAntiDeuteron_Event_nParticles,"AntiDeuteron_Event_nParticles/s");




  PostData(1,fSaveTree_Lambda);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiLambda);
  PostData(4,fSaveTree_AntiDeuteron);
  PostData(5,fHistoList);



} // end of UserCreateOutputObjects









void AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec(Option_t*)
{

//  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAODEvent)::Fatal("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AOD event found!");

  fHeader = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  if(!fHeader)::Fatal("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No Header found!");

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse)::Fatal("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No PIDResponse found!");

  if(fIsMC == true){
  
    fMCEvent = MCEvent(); 
    if(!fMCEvent)::Fatal("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No MC event found!");

  }


  // debug analysis
  bool DebugEventSelection  = false;
  bool DebugMC = false;

  // define event cuts
  double PrimaryVertexMaxZ  = 10.0; // cm
  double Centrality_min	    = 0.0;
  double Centrality_max	    = 100.0;

  if(fCollisionSystem == 1) // Central collisions
  {
    Centrality_min = 0.0;
    Centrality_max = 10.0;
  }

  if(fCollisionSystem == 2) // Semi-Central collisions
  {
    Centrality_min = 30.0;
    Centrality_max = 50.0;
  }

  // use only events containing tracks
  int nTracks = fAODEvent->GetNumberOfTracks();
  if(nTracks == 0) return;

  // use only events containing v0s
  int nV0s = fAODEvent->GetNumberOfV0s();
  if(nV0s == 0) return;

  // get primary vertex
  AliAODVertex *PrimaryVertex = fAODEvent->GetPrimaryVertex();
  if(!PrimaryVertex)::Warning("AliAnalsisTask_Ld_CreateTrees_PairsOnlyd::UserExec","No AliAODVertex object found!");
  double PrimaryVertexPos[3] = {-999.0,-999.0,-999.0};
  PrimaryVertex->GetXYZ(PrimaryVertexPos);

  // apply cut on z-position of primary vertex
  double PrimaryVertexZ = PrimaryVertexPos[2]; // cm
  if(TMath::IsNaN(PrimaryVertexZ)) return;
  if(TMath::Abs(PrimaryVertexZ) > PrimaryVertexMaxZ) return;

  // apply centrality cut
  double Centrality = -999.0;
  AliMultSelection *MultSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
  Centrality = MultSelection->GetMultiplicityPercentile("V0M");
  if(TMath::IsNaN(Centrality)) return;

  if((fCollisionSystem == 1) || (fCollisionSystem == 2)){
    if((Centrality < Centrality_min) || (Centrality > Centrality_max)) return;
  }


  //combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8
  unsigned int Multiplicity = fHeader->GetRefMultiplicityComb08();



  // get event information
  unsigned int PeriodNumber	  = fAODEvent->GetPeriodNumber();
  unsigned int OrbitNumber	  = fAODEvent->GetOrbitNumber();
  unsigned short BunchCrossNumber = fAODEvent->GetBunchCrossNumber();
  int RunNumber			  = fAODEvent->GetRunNumber();
  float BField			  = fAODEvent->GetMagneticField();
  if(TMath::IsNaN(BField)) return;
  unsigned long EventID = 0;
  unsigned long Seed = 0;
  int RandomNumber = 0;
  int nTracksMC = 0;


  // EventID -> https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventInfo
  if(fIsMC == false) EventID = (unsigned long)BunchCrossNumber + ((unsigned long)OrbitNumber*3564) + ((unsigned long)PeriodNumber*16777215*3564);

  // EventID (MC)
  if(fIsMC == true){

    TRandom3 *RandomGenerator = new TRandom3();
    RandomGenerator->SetSeed(0);
    Seed = RandomGenerator->GetSeed();
    unsigned int Offset1 = 10000000;
    unsigned int Offset2 = 100000;
    RandomNumber = RandomGenerator->Integer(Offset2);

    int IntegerPrimaryVertexZ = TMath::Abs((int)(std::trunc(PrimaryVertexZ*1000000)));
    nTracksMC = fMCEvent->GetNumberOfTracks();
    EventID = ((((unsigned long) RunNumber * Offset1) + (unsigned long)IntegerPrimaryVertexZ) * Offset2 * 10) + (unsigned long)RandomNumber;

  } // end of fIsMC == true

                                                                                                                                               
  TF1 *FitPionDCAxy_min = new TF1("FitPionDCAxy_min","pol2",0.0,6.0);
  FitPionDCAxy_min->FixParameter(0,0.0635603);
  FitPionDCAxy_min->FixParameter(1,-0.0593139);
  FitPionDCAxy_min->FixParameter(2,0.0171978);
  
  TF1 *FitPionDCAz_min = new TF1("FitPionDCAz_min","pol2",0.0,6.0);
  FitPionDCAz_min->FixParameter(0,0.071036);
  FitPionDCAz_min->FixParameter(1,-0.0441316);
  FitPionDCAz_min->FixParameter(2,0.00701049);





  // print event information
  if(DebugEventSelection)
  {

    cout << "" << endl;
    cout << "fCollisionSystem:\t\t" << fCollisionSystem << std::endl;
    cout << "PeriodNumber:\t\t\t" << PeriodNumber << endl;
    cout << "RunNumber:\t\t\t" << RunNumber << endl;
    cout << "OrbitNumber:\t\t\t" << OrbitNumber << endl;
    cout << "BunchCrossNumber:\t\t" << BunchCrossNumber << endl;
    cout << "Unique Event ID:\t\t" << EventID << endl;
    cout << "Centrality:\t" << Centrality << " %" << endl;
    cout << "Multiplicity:\t" << Multiplicity << endl;
    cout << "Number of tracks in event:\t" << nTracks << endl;
    cout << "Number of tracks in MC event:\t" << nTracksMC << endl;
    cout << "Seed of RandomGenerator:\t" << Seed << endl;
    cout << "RandomNumber:\t\t\t" << RandomNumber << endl;
    cout << "z-position of primary vertex:\t" << PrimaryVertexZ << " cm" << endl;

  } // end of DebugEventSelection



  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ proton selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    float Lambda_px;
    float Lambda_py;
    float Lambda_pz;
    float Lambda_px_Generated = 0.0;
    float Lambda_py_Generated = 0.0;
    float Lambda_pz_Generated = 0.0;
    float Lambda_Eta;
    float Lambda_Phi;
    float Lambda_TransverseRadius;
    float Lambda_CosinePointingAngle;
    float Lambda_DCAv0ToPrimaryVertex;
    float Lambda_DCAv0Daughters;
    float Lambda_Alpha;
    float Lambda_qT;
    float Lambda_DecayLength;
    float Lambda_OpenAngle;
    int	  Lambda_PDG_Daughter1;
    int	  Lambda_PDG_Daughter2;
    int	  Lambda_PDG_v01;
    int	  Lambda_PDG_v02;
    int	  Lambda_PDG_Mother1;
    int	  Lambda_PDG_Mother2;
    bool  Lambda_SameV0;

  float     Lambda_Daughter_Proton_px;
  float     Lambda_Daughter_Proton_py;
  float     Lambda_Daughter_Proton_pz;
  float     Lambda_Daughter_Proton_px_Generated;
  float     Lambda_Daughter_Proton_py_Generated;
  float     Lambda_Daughter_Proton_pz_Generated;
  float     Lambda_Daughter_Proton_px_DecayVertex;
  float     Lambda_Daughter_Proton_py_DecayVertex;
  float     Lambda_Daughter_Proton_pz_DecayVertex;
  float     Lambda_Daughter_Proton_pTPC;
  float     Lambda_Daughter_Proton_Eta;
  float     Lambda_Daughter_Proton_Phi;
  float     Lambda_Daughter_Proton_TPC_Chi2;
  float     Lambda_Daughter_Proton_TPC_dEdx;
  float     Lambda_Daughter_Proton_TPC_dEdx_nSigma;
  float     Lambda_Daughter_Proton_TOF_Mass2;
  float     Lambda_Daughter_Proton_TOF_Mass2_nSigma;
  float     Lambda_Daughter_Proton_ITS_dEdx;
  float     Lambda_Daughter_Proton_ITS_dEdx_nSigma;
  float     Lambda_Daughter_Proton_DCAxy;
  float     Lambda_Daughter_Proton_DCAz;
  unsigned short    Lambda_Daughter_Proton_TPC_nCrossedRows;
  unsigned short    Lambda_Daughter_Proton_TPC_nSharedCluster;
  unsigned short    Lambda_Daughter_Proton_TPC_nFindableCluster;
  unsigned short    Lambda_Daughter_Proton_TPC_nCluster;
  unsigned short    Lambda_Daughter_Proton_ITS_nCluster;
  unsigned int      Lambda_Daughter_Proton_ID = 0;

  float     Lambda_Daughter_AntiPion_px;
  float     Lambda_Daughter_AntiPion_py;
  float     Lambda_Daughter_AntiPion_pz;
  float     Lambda_Daughter_AntiPion_px_Generated;
  float     Lambda_Daughter_AntiPion_py_Generated;
  float     Lambda_Daughter_AntiPion_pz_Generated;
  float     Lambda_Daughter_AntiPion_px_DecayVertex;
  float     Lambda_Daughter_AntiPion_py_DecayVertex;
  float     Lambda_Daughter_AntiPion_pz_DecayVertex;
  float     Lambda_Daughter_AntiPion_pTPC;
  float     Lambda_Daughter_AntiPion_Eta;
  float     Lambda_Daughter_AntiPion_Phi;
  float     Lambda_Daughter_AntiPion_TPC_Chi2;
  float     Lambda_Daughter_AntiPion_TPC_dEdx;
  float     Lambda_Daughter_AntiPion_TPC_dEdx_nSigma;
  float     Lambda_Daughter_AntiPion_TOF_Mass2;
  float     Lambda_Daughter_AntiPion_TOF_Mass2_nSigma;
  float     Lambda_Daughter_AntiPion_ITS_dEdx;
  float     Lambda_Daughter_AntiPion_ITS_dEdx_nSigma;
  float     Lambda_Daughter_AntiPion_DCAxy;
  float     Lambda_Daughter_AntiPion_DCAz;
  unsigned short    Lambda_Daughter_AntiPion_TPC_nCrossedRows;
  unsigned short    Lambda_Daughter_AntiPion_TPC_nSharedCluster;
  unsigned short    Lambda_Daughter_AntiPion_TPC_nFindableCluster;
  unsigned short    Lambda_Daughter_AntiPion_TPC_nCluster;
  unsigned short    Lambda_Daughter_AntiPion_ITS_nCluster;
  unsigned int      Lambda_Daughter_AntiPion_ID = 0;


  TTree *fTempTree_Lambda = new TTree("fTempTree_Lambda","fTempTree_Lambda");
  fTempTree_Lambda->Branch("Lambda_px",&Lambda_px,"Lambda_px/F");
  fTempTree_Lambda->Branch("Lambda_py",&Lambda_py,"Lambda_py/F");
  fTempTree_Lambda->Branch("Lambda_pz",&Lambda_pz,"Lambda_pz/F");
  fTempTree_Lambda->Branch("Lambda_px_Generated",&Lambda_px_Generated,"Lambda_px_Generated/F");
  fTempTree_Lambda->Branch("Lambda_py_Generated",&Lambda_py_Generated,"Lambda_py_Generated/F");
  fTempTree_Lambda->Branch("Lambda_pz_Generated",&Lambda_pz_Generated,"Lambda_pz_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Eta",&Lambda_Eta,"Lambda_Eta/F");
  fTempTree_Lambda->Branch("Lambda_Phi",&Lambda_Phi,"Lambda_Phi/F");
  fTempTree_Lambda->Branch("Lambda_TransverseRadius",&Lambda_TransverseRadius,"Lambda_TransverseRadius/F");
  fTempTree_Lambda->Branch("Lambda_CosinePointingAngle",&Lambda_CosinePointingAngle,"Lambda_CosinePointingAngle/F");
  fTempTree_Lambda->Branch("Lambda_DCAv0ToPrimaryVertex",&Lambda_DCAv0ToPrimaryVertex,"Lambda_DCAv0ToPrimaryVertex/F");
  fTempTree_Lambda->Branch("Lambda_DCAv0Daughters",&Lambda_DCAv0Daughters,"Lambda_DCAv0Daughters/F");
  fTempTree_Lambda->Branch("Lambda_Alpha",&Lambda_Alpha,"Lambda_Alpha/F");
  fTempTree_Lambda->Branch("Lambda_qT",&Lambda_qT,"Lambda_qT/F");
  fTempTree_Lambda->Branch("Lambda_DecayLength",&Lambda_DecayLength,"Lambda_DecayLength/F");
  fTempTree_Lambda->Branch("Lambda_OpenAngle",&Lambda_OpenAngle,"Lambda_OpenAngle/F");
  fTempTree_Lambda->Branch("Lambda_PDG_Daughter1",&Lambda_PDG_Daughter1,"Lambda_PDG_Daughter1/I");
  fTempTree_Lambda->Branch("Lambda_PDG_Daughter2",&Lambda_PDG_Daughter2,"Lambda_PDG_Daughter2/I");
  fTempTree_Lambda->Branch("Lambda_PDG_v01",&Lambda_PDG_v01,"Lambda_PDG_v01/I");
  fTempTree_Lambda->Branch("Lambda_PDG_v02",&Lambda_PDG_v02,"Lambda_PDG_v02/I");
  fTempTree_Lambda->Branch("Lambda_PDG_Mother1",&Lambda_PDG_Mother1,"Lambda_PDG_Mother1/I");
  fTempTree_Lambda->Branch("Lambda_PDG_Mother2",&Lambda_PDG_Mother2,"Lambda_PDG_Mother2/I");
  fTempTree_Lambda->Branch("Lambda_SameV0",&Lambda_SameV0,"Lambda_SameV0/O");

  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_px",&Lambda_Daughter_Proton_px,"Lambda_Daughter_Proton_px/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_py",&Lambda_Daughter_Proton_py,"Lambda_Daughter_Proton_py/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_pz",&Lambda_Daughter_Proton_pz,"Lambda_Daughter_Proton_pz/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_px_Generated",&Lambda_Daughter_Proton_px_Generated,"Lambda_Daughter_Proton_px_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_py_Generated",&Lambda_Daughter_Proton_py_Generated,"Lambda_Daughter_Proton_py_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_pz_Generated",&Lambda_Daughter_Proton_pz_Generated,"Lambda_Daughter_Proton_pz_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_px_DecayVertex",&Lambda_Daughter_Proton_px_DecayVertex,"Lambda_Daughter_Proton_px_DecayVertex/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_py_DecayVertex",&Lambda_Daughter_Proton_py_DecayVertex,"Lambda_Daughter_Proton_py_DecayVertex/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_pz_DecayVertex",&Lambda_Daughter_Proton_pz_DecayVertex,"Lambda_Daughter_Proton_pz_DecayVertex/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_pTPC",&Lambda_Daughter_Proton_pTPC,"Lambda_Daughter_Proton_pTPC/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_Eta",&Lambda_Daughter_Proton_Eta,"Lambda_Daughter_Proton_Eta/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_Phi",&Lambda_Daughter_Proton_Phi,"Lambda_Daughter_Proton_Phi/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_Chi2",&Lambda_Daughter_Proton_TPC_Chi2,"Lambda_Daughter_Proton_TPC_Chi2/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_dEdx",&Lambda_Daughter_Proton_TPC_dEdx,"Lambda_Daughter_Proton_TPC_dEdx/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_dEdx_nSigma",&Lambda_Daughter_Proton_TPC_dEdx_nSigma,"Lambda_Daughter_Proton_TPC_dEdx_nSigma/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TOF_Mass2",&Lambda_Daughter_Proton_TOF_Mass2,"Lambda_Daughter_Proton_TOF_Mass2/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TOF_Mass2_nSigma",&Lambda_Daughter_Proton_TOF_Mass2_nSigma,"Lambda_Daughter_Proton_TOF_Mass2_nSigma/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_dEdx",&Lambda_Daughter_Proton_ITS_dEdx,"Lambda_Daughter_Proton_ITS_dEdx/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_dEdx_nSigma",&Lambda_Daughter_Proton_ITS_dEdx_nSigma,"Lambda_Daughter_Proton_ITS_dEdx_nSigma/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_DCAxy",&Lambda_Daughter_Proton_DCAxy,"Lambda_Daughter_Proton_DCAxy/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_DCAz",&Lambda_Daughter_Proton_DCAz,"Lambda_Daughter_Proton_DCAz/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nCrossedRows",&Lambda_Daughter_Proton_TPC_nCrossedRows,"Lambda_Daughter_Proton_TPC_nCrossedRows/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nSharedCluster",&Lambda_Daughter_Proton_TPC_nSharedCluster,"Lambda_Daughter_Proton_TPC_nSharedCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nFindableCluster",&Lambda_Daughter_Proton_TPC_nFindableCluster,"Lambda_Daughter_Proton_TPC_nFindableCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nCluster",&Lambda_Daughter_Proton_TPC_nCluster,"Lambda_Daughter_Proton_TPC_nCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_nCluster",&Lambda_Daughter_Proton_ITS_nCluster,"Lambda_Daughter_Proton_ITS_nCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_Proton_ID",&Lambda_Daughter_Proton_ID,"Lambda_Daughter_Proton_ID/i");

  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_px",&Lambda_Daughter_AntiPion_px,"Lambda_Daughter_AntiPion_px/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_py",&Lambda_Daughter_AntiPion_py,"Lambda_Daughter_AntiPion_py/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz",&Lambda_Daughter_AntiPion_pz,"Lambda_Daughter_AntiPion_pz/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_px_Generated",&Lambda_Daughter_AntiPion_px_Generated,"Lambda_Daughter_AntiPion_px_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_py_Generated",&Lambda_Daughter_AntiPion_py_Generated,"Lambda_Daughter_AntiPion_py_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz_Generated",&Lambda_Daughter_AntiPion_pz_Generated,"Lambda_Daughter_AntiPion_pz_Generated/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_px_DecayVertex",&Lambda_Daughter_AntiPion_px_DecayVertex,"Lambda_Daughter_AntiPion_px_DecayVertex/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_py_DecayVertex",&Lambda_Daughter_AntiPion_py_DecayVertex,"Lambda_Daughter_AntiPion_py_DecayVertex/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz_DecayVertex",&Lambda_Daughter_AntiPion_pz_DecayVertex,"Lambda_Daughter_AntiPion_pz_DecayVertex/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_pTPC",&Lambda_Daughter_AntiPion_pTPC,"Lambda_Daughter_AntiPion_pTPC/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_Eta",&Lambda_Daughter_AntiPion_Eta,"Lambda_Daughter_AntiPion_Eta/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_Phi",&Lambda_Daughter_AntiPion_Phi,"Lambda_Daughter_AntiPion_Phi/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_Chi2",&Lambda_Daughter_AntiPion_TPC_Chi2,"Lambda_Daughter_AntiPion_TPC_Chi2/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_dEdx",&Lambda_Daughter_AntiPion_TPC_dEdx,"Lambda_Daughter_AntiPion_TPC_dEdx/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_dEdx_nSigma",&Lambda_Daughter_AntiPion_TPC_dEdx_nSigma,"Lambda_Daughter_AntiPion_TPC_dEdx_nSigma/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TOF_Mass2",&Lambda_Daughter_AntiPion_TOF_Mass2,"Lambda_Daughter_AntiPion_TOF_Mass2/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TOF_Mass2_nSigma",&Lambda_Daughter_AntiPion_TOF_Mass2_nSigma,"Lambda_Daughter_AntiPion_TOF_Mass2_nSigma/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_dEdx",&Lambda_Daughter_AntiPion_ITS_dEdx,"Lambda_Daughter_AntiPion_ITS_dEdx/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_dEdx_nSigma",&Lambda_Daughter_AntiPion_ITS_dEdx_nSigma,"Lambda_Daughter_AntiPion_ITS_dEdx_nSigma/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_DCAxy",&Lambda_Daughter_AntiPion_DCAxy,"Lambda_Daughter_AntiPion_DCAxy/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_DCAz",&Lambda_Daughter_AntiPion_DCAz,"Lambda_Daughter_AntiPion_DCAz/F");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nCrossedRows",&Lambda_Daughter_AntiPion_TPC_nCrossedRows,"Lambda_Daughter_AntiPion_TPC_nCrossedRows/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nSharedCluster",&Lambda_Daughter_AntiPion_TPC_nSharedCluster,"Lambda_Daughter_AntiPion_TPC_nSharedCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nFindableCluster",&Lambda_Daughter_AntiPion_TPC_nFindableCluster,"Lambda_Daughter_AntiPion_TPC_nFindableCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nCluster",&Lambda_Daughter_AntiPion_TPC_nCluster,"Lambda_Daughter_AntiPion_TPC_nCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_nCluster",&Lambda_Daughter_AntiPion_ITS_nCluster,"Lambda_Daughter_AntiPion_ITS_nCluster/s");
  fTempTree_Lambda->Branch("Lambda_Daughter_AntiPion_ID",&Lambda_Daughter_AntiPion_ID,"Lambda_Daughter_AntiPion_ID/i");




  unsigned short nLambdasSelected = 0;

  for(int V0 = 0; V0 < nV0s; V0++)  // Lambda loop
  { 

    AliAODv0 *v0 = dynamic_cast<AliAODv0*>(fAODEvent->GetV0(V0));
    if(!v0)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODv0 found");
    if(!v0) continue;

    AliAODTrack *ProtonTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    if(!ProtonTrack) continue;

    AliAODTrack *AntiPionTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    if(!AntiPionTrack) continue;

    bool PassedProtonCuts = CheckProtonCuts(*ProtonTrack,*fPIDResponse,true,RunNumber);
    if(!PassedProtonCuts) continue;

    double PionDCAxy_min = FitPionDCAxy_min->Eval(AntiPionTrack->Pt());
    double PionDCAz_min = FitPionDCAz_min->Eval(AntiPionTrack->Pt());
    bool PassedPionCuts = CheckPionCuts(*AntiPionTrack,*fPIDResponse,false,RunNumber,PionDCAxy_min,PionDCAz_min);
    if(!PassedPionCuts) continue;


    AliAODVertex *DecayVertex = v0->GetSecondaryVtx();

    AliExternalTrackParam *ProtonTrackParam = new AliExternalTrackParam();
    ProtonTrackParam->CopyFromVTrack(ProtonTrack);
    double Proton_DCA_new[2];
    double Proton_CovarianceMatrix_new[3];
    ProtonTrackParam->PropagateToDCA(DecayVertex,BField,10.0,Proton_DCA_new,Proton_CovarianceMatrix_new);
    double ProtonMomentumPropagated[3];
    ProtonTrackParam->GetPxPyPz(ProtonMomentumPropagated);

    AliExternalTrackParam *AntiPionTrackParam = new AliExternalTrackParam();
    AntiPionTrackParam->CopyFromVTrack(AntiPionTrack);
    double AntiPion_DCA_new[2];
    double AntiPion_CovarianceMatrix_new[3];
    AntiPionTrackParam->PropagateToDCA(DecayVertex,BField,10.0,AntiPion_DCA_new,AntiPion_CovarianceMatrix_new);
    double AntiPionMomentumPropagated[3];
    AntiPionTrackParam->GetPxPyPz(AntiPionMomentumPropagated);


    bool PrintDecayMomentaOnScreen = false;
    if(PrintDecayMomentaOnScreen == true){

      std::cout << "\nProtonTrack->Px() = " << ProtonTrack->Px() << "\t ProtonTrack->Py() = " << ProtonTrack->Py() << "\t ProtonTrack->Pz() = " << ProtonTrack->Pz() << std::endl;
      std::cout << "v0->MomPosX() = " << v0->MomPosX() << "\t v0->MomPosY() = " << v0->MomPosY() << "\tv0->MomPosZ() = " << v0->MomPosZ() << std::endl;
      std::cout << "px_prop = " << ProtonMomentumPropagated[0] << "\t\t py_prop = " << ProtonMomentumPropagated[1] << "\t\t pz_prop = " << ProtonMomentumPropagated[2] << std::endl;
    std::cout << "AntiPionTrack->Px() = " << AntiPionTrack->Px() << "\t AntiPionTrack->Py() = " << AntiPionTrack->Py() << "\tAntiPionTrack->Pz() = " << AntiPionTrack->Pz() << std::endl;
    std::cout << "v0->MomNegX() = " << v0->MomNegX() << "\t v0->MomNegY() = " << v0->MomNegY() << "\tv0->MomNegZ() = " << v0->MomNegZ() << std::endl;
    std::cout << "px_prop = " << AntiPionMomentumPropagated[0] << "\t py_prop = " << AntiPionMomentumPropagated[1] << "\t pz_prop = " << AntiPionMomentumPropagated[2] << std::endl;

    } // end of PrintDecayMomentaOnScreen


    double MomentumDaughterPosX = ProtonMomentumPropagated[0];
    double MomentumDaughterPosY = ProtonMomentumPropagated[1];
    double MomentumDaughterPosZ = ProtonMomentumPropagated[2];
    double MomentumDaughterNegX = AntiPionMomentumPropagated[0];
    double MomentumDaughterNegY = AntiPionMomentumPropagated[1];
    double MomentumDaughterNegZ = AntiPionMomentumPropagated[2];
    double MomentumV0X = v0->MomV0X();
    double MomentumV0Y = v0->MomV0Y();
    double MomentumV0Z = v0->MomV0Z();


    // apply Lambda cuts
    double MassInvariant = CalculateInvariantMassLambda(ProtonMomentumPropagated,AntiPionMomentumPropagated,true);
    double WrongMassInvariant = CalculateInvariantMassLambda(AntiPionMomentumPropagated,ProtonMomentumPropagated,true);
    bool PassedLambdaCuts = CheckLambdaCuts(*v0,PrimaryVertexPos,*fPIDResponse,true,RunNumber,MassInvariant,WrongMassInvariant);
    if(!PassedLambdaCuts) continue;




    float Generated_px = -999.0;
    float Generated_py = -999.0;
    float Generated_pz = -999.0;
    float Generated_px_Daughter1 = -999.0;
    float Generated_py_Daughter1 = -999.0;
    float Generated_pz_Daughter1 = -999.0;
    float Generated_px_Daughter2 = -999.0;
    float Generated_py_Daughter2 = -999.0;
    float Generated_pz_Daughter2 = -999.0;
    int PDG_Daughter1 = 0;
    int PDG_Daughter2 = 0;
    int PDG_v01 = 0;
    int PDG_v02 = 0;
    int PDG_Mother1 = 0;
    int PDG_Mother2 = 0;
    bool SameV0 = false;

    int Label_Daughter1 = TMath::Abs(ProtonTrack->GetLabel());
    int Label_Daughter2 = TMath::Abs(AntiPionTrack->GetLabel());


    if(fIsMC == true){

      // get list of generated particles
      TClonesArray *MCArray = dynamic_cast <TClonesArray*> (fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!MCArray) continue;

      // Daughter1, v01, Mother1
      AliMCParticle *MCParticleDaughter1 = (AliMCParticle*) MCArray->At(Label_Daughter1);
      PDG_Daughter1 = MCParticleDaughter1->PdgCode();

      int Label_v01 = TMath::Abs(MCParticleDaughter1->GetMother());
      AliMCParticle *MCv01 = (AliMCParticle*) MCArray->At(Label_v01);
      PDG_v01 = MCv01->PdgCode();

      int Label_Mother1 = TMath::Abs(MCv01->GetMother());
      AliMCParticle *MCParticleMother1 = (AliMCParticle*) MCArray->At(Label_Mother1);
      PDG_Mother1 = MCParticleMother1->PdgCode();


      // Daughter2, v02, Mother2
      AliMCParticle *MCParticleDaughter2 = (AliMCParticle*) fMCEvent->GetTrack(Label_Daughter2);
      PDG_Daughter2 = MCParticleDaughter2->PdgCode();

      int Label_v02 = TMath::Abs(MCParticleDaughter2->GetMother());
      AliMCParticle *MCv02 = (AliMCParticle*) MCArray->At(Label_v02);
      PDG_v02 = MCv02->PdgCode();

      int Label_Mother2 = TMath::Abs(MCv02->GetMother());
      AliMCParticle *MCParticleMother2 = (AliMCParticle*) MCArray->At(Label_Mother2);
      PDG_Mother2 = MCParticleMother2->PdgCode();



      // check if daughters come from the same v0 
      SameV0 = false;
      if(Label_v01 == Label_v02) SameV0 = true;


      // where do the particles come from?
      if(MCParticleDaughter1->IsPhysicalPrimary() == true)	  PDG_v01 = 1;
      if(MCParticleDaughter1->IsSecondaryFromMaterial() == true)  PDG_v01 = 2;

      if(MCParticleDaughter2->IsPhysicalPrimary() == true)	  PDG_v02 = 1;
      if(MCParticleDaughter2->IsSecondaryFromMaterial() == true)  PDG_v02 = 2;

      if(MCv01->IsPhysicalPrimary() == true)	    PDG_Mother1 = 1;
      if(MCv01->IsSecondaryFromMaterial() == true)  PDG_Mother1 = 2;

      if(MCv02->IsPhysicalPrimary() == true)	    PDG_Mother2 = 1;
      if(MCv02->IsSecondaryFromMaterial() == true)  PDG_Mother2 = 2;

      if(SameV0 == true){

	// save generated momenta
	Generated_px = MCv01->Px();
	Generated_py = MCv01->Py();
	Generated_pz = MCv01->Pz();

      }

      Generated_px_Daughter1 = MCParticleDaughter1->Px();
      Generated_py_Daughter1 = MCParticleDaughter1->Py();
      Generated_pz_Daughter1 = MCParticleDaughter1->Pz();

      Generated_px_Daughter2 = MCParticleDaughter2->Px();
      Generated_py_Daughter2 = MCParticleDaughter2->Py();
      Generated_pz_Daughter2 = MCParticleDaughter2->Pz();


      if(DebugMC == true){

	std::cout << "\n v0: " << V0 << "/" << nV0s << std::endl;
	std::cout << "Label_Daughter1: " << Label_Daughter1 << "\t PDG_Daughter1: " << PDG_Daughter1 << "\t Label_Daughter2: " << Label_Daughter2 << "\t PDG_Daughter2: " << PDG_Daughter2 << std::endl;
	std::cout << "Label_v01: " << Label_v01 << "\t PDG_v01: " << PDG_v01 << "\t Label_v02: " << Label_v02 << "\t PDG_v02: " << PDG_v02 << "\t SameV0: " << SameV0 << std::endl;
	std::cout << "Label_Mother1: " << Label_Mother1 << "\t PDG_Mother1: " << PDG_Mother1 << "\t Label_Mother2: " << Label_Mother2 << "\t PDG_Mother2: " << PDG_Mother2 << "\n" << std::endl;

      }


    } // end of fIsMC == true






  
    TVector3 *MomentumDaughterPositive = new TVector3();
    TVector3 *MomentumDaughterNegative = new TVector3();
    TVector3 *MomentumV0 = new TVector3();
    MomentumDaughterPositive->SetXYZ(MomentumDaughterPosX,MomentumDaughterPosY,MomentumDaughterPosZ);
    MomentumDaughterNegative->SetXYZ(MomentumDaughterNegX,MomentumDaughterNegY,MomentumDaughterNegZ);
    MomentumV0->SetXYZ(MomentumV0X,MomentumV0Y,MomentumV0Z);
  
    double p_Longitudinal_Positive = MomentumDaughterPositive->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of positively charged daughter
    double p_Longitudinal_Negative = MomentumDaughterNegative->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of negatively charged daughter
  
    double Alpha = (p_Longitudinal_Positive - p_Longitudinal_Negative) / (p_Longitudinal_Positive + p_Longitudinal_Negative);
    double qT = MomentumDaughterNegative->Perp(*MomentumV0);

    float Proton_xv[2];
    float Proton_yv[3];
    ProtonTrack->GetImpactParameters(Proton_xv,Proton_yv);
    float Proton_DCAxy = Proton_xv[0];
    float Proton_DCAz = Proton_xv[1];

    double Proton_TOF_m2	   = -999.0;
    double Proton_TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus ProtonStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,ProtonTrack);
    bool ProtonTOFisOK = false;
    if(ProtonStatusTOF == AliPIDResponse::kDetPidOk) ProtonTOFisOK = true;
    
    if(ProtonTOFisOK == true){

      Proton_TOF_m2	    = CalculateMassSquareTOF(*ProtonTrack);
      Proton_TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(ProtonTrack->Pt(),Proton_TOF_m2,1,RunNumber);

    }

    double Proton_ITS_dEdx	    = -999.0;
    double Proton_ITS_dEdx_nSigma   = -999.0;
    int Proton_ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus ProtonStatusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,ProtonTrack);
    bool ProtonITSisOK = false;
    if(ProtonStatusITS == AliPIDResponse::kDetPidOk) ProtonITSisOK = true;
    
    if(ProtonITSisOK == true){

      Proton_ITS_dEdx	      = ProtonTrack->GetITSsignal();
      Proton_ITS_dEdx_nSigma  = CalculateSigmadEdxITS(*ProtonTrack,1,RunNumber);
      Proton_ITS_nCluster     = ProtonTrack->GetITSNcls();
      if(Proton_ITS_nCluster < 0) Proton_ITS_nCluster = 0;

    }



    float AntiPion_xv[2];
    float AntiPion_yv[3];
    AntiPionTrack->GetImpactParameters(AntiPion_xv,AntiPion_yv);
    float AntiPion_DCAxy = AntiPion_xv[0];
    float AntiPion_DCAz = AntiPion_xv[1];

    double AntiPion_TOF_m2	   = -999.0;
    double AntiPion_TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus AntiPionStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,AntiPionTrack);
    bool AntiPionTOFisOK = false;
    if(AntiPionStatusTOF == AliPIDResponse::kDetPidOk) AntiPionTOFisOK = true;
    
    if(AntiPionTOFisOK == true){

      AntiPion_TOF_m2	      = CalculateMassSquareTOF(*AntiPionTrack);
      AntiPion_TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(AntiPionTrack->Pt(),AntiPion_TOF_m2,6,RunNumber);

    }

    double AntiPion_ITS_dEdx	    = -999.0;
    double AntiPion_ITS_dEdx_nSigma   = -999.0;
    int AntiPion_ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus AntiPionStatusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,AntiPionTrack);
    bool AntiPionITSisOK = false;
    if(AntiPionStatusITS == AliPIDResponse::kDetPidOk) AntiPionITSisOK = true;
    
    if(AntiPionITSisOK == true){

      AntiPion_ITS_dEdx	      = AntiPionTrack->GetITSsignal();
      AntiPion_ITS_dEdx_nSigma  = CalculateSigmadEdxITS(*AntiPionTrack,6,RunNumber);
      AntiPion_ITS_nCluster     = AntiPionTrack->GetITSNcls();
      if(AntiPion_ITS_nCluster < 0) AntiPion_ITS_nCluster = 0;

    }




    Lambda_px			  = (float)v0->MomV0X();
    Lambda_py			  = (float)v0->MomV0Y();
    Lambda_pz			  = (float)v0->MomV0Z();
    Lambda_px_Generated		  = Generated_px;
    Lambda_py_Generated		  = Generated_py;
    Lambda_pz_Generated		  = Generated_pz;
    Lambda_Eta			  = (float)v0->PseudoRapV0();
    Lambda_Phi			  = (float)v0->Phi();
    Lambda_TransverseRadius	  = (float)v0->RadiusV0();
    Lambda_CosinePointingAngle	  = (float)TMath::Abs(v0->CosPointingAngle(PrimaryVertexPos));
    Lambda_DCAv0ToPrimaryVertex	  = (float)v0->DcaV0ToPrimVertex();
    Lambda_DCAv0Daughters	  = (float)v0->DcaV0Daughters();
    Lambda_Alpha		  = (float)Alpha;
    Lambda_qT			  = (float)qT;
    Lambda_DecayLength		  = (float)v0->DecayLengthV0(PrimaryVertexPos);
    Lambda_OpenAngle		  = (float)v0->OpenAngleV0();
    Lambda_PDG_Daughter1	  = PDG_Daughter1;
    Lambda_PDG_Daughter2	  = PDG_Daughter2;
    Lambda_PDG_v01		  = PDG_v01;
    Lambda_PDG_v02		  = PDG_v02;
    Lambda_PDG_Mother1		  = PDG_Mother1;
    Lambda_PDG_Mother2		  = PDG_Mother2;
    Lambda_SameV0		  = SameV0;

    Lambda_Daughter_Proton_px			    = ProtonTrack->Px();
    Lambda_Daughter_Proton_py			    = ProtonTrack->Py();
    Lambda_Daughter_Proton_pz			    = ProtonTrack->Pz();
    Lambda_Daughter_Proton_px_Generated		    = Generated_px_Daughter1;
    Lambda_Daughter_Proton_py_Generated		    = Generated_py_Daughter1;
    Lambda_Daughter_Proton_pz_Generated		    = Generated_pz_Daughter1;
    Lambda_Daughter_Proton_px_DecayVertex	    = MomentumDaughterPosX;
    Lambda_Daughter_Proton_py_DecayVertex	    = MomentumDaughterPosY;
    Lambda_Daughter_Proton_pz_DecayVertex	    = MomentumDaughterPosZ;
    Lambda_Daughter_Proton_pTPC			    = ProtonTrack->GetTPCmomentum();
    Lambda_Daughter_Proton_Eta			    = ProtonTrack->Eta();
    Lambda_Daughter_Proton_Phi			    = ProtonTrack->Phi();
    Lambda_Daughter_Proton_TPC_Chi2		    = ProtonTrack->GetTPCchi2();
    Lambda_Daughter_Proton_TPC_dEdx		    = ProtonTrack->GetTPCsignal();
    Lambda_Daughter_Proton_TPC_dEdx_nSigma	    = (float)fPIDResponse->NumberOfSigmasTPC(ProtonTrack,AliPID::kProton);
    Lambda_Daughter_Proton_TOF_Mass2		    = (float)Proton_TOF_m2;
    Lambda_Daughter_Proton_TOF_Mass2_nSigma	    = (float)Proton_TOF_m2_nSigma;
    Lambda_Daughter_Proton_ITS_dEdx		    = (float)Proton_ITS_dEdx;
    Lambda_Daughter_Proton_ITS_dEdx_nSigma	    = (float)Proton_ITS_dEdx_nSigma;
    Lambda_Daughter_Proton_DCAxy		    = Proton_DCAxy;
    Lambda_Daughter_Proton_DCAz			    = Proton_DCAz;
    Lambda_Daughter_Proton_TPC_nCrossedRows	    = ProtonTrack->GetTPCCrossedRows();
    Lambda_Daughter_Proton_TPC_nSharedCluster	    = ProtonTrack->GetTPCnclsS();
    Lambda_Daughter_Proton_TPC_nFindableCluster	    = ProtonTrack->GetTPCNclsF();
    Lambda_Daughter_Proton_TPC_nCluster		    = ProtonTrack->GetTPCNcls();
    Lambda_Daughter_Proton_ITS_nCluster		    = (unsigned short)Proton_ITS_nCluster;
    Lambda_Daughter_Proton_ID			    = Label_Daughter1;

    Lambda_Daughter_AntiPion_px			    = AntiPionTrack->Px();
    Lambda_Daughter_AntiPion_py			    = AntiPionTrack->Py();
    Lambda_Daughter_AntiPion_pz			    = AntiPionTrack->Pz();
    Lambda_Daughter_AntiPion_px_Generated	    = Generated_px_Daughter2;
    Lambda_Daughter_AntiPion_py_Generated	    = Generated_py_Daughter2;
    Lambda_Daughter_AntiPion_pz_Generated	    = Generated_pz_Daughter2;
    Lambda_Daughter_AntiPion_px_DecayVertex	    = MomentumDaughterNegX;
    Lambda_Daughter_AntiPion_py_DecayVertex	    = MomentumDaughterNegY;
    Lambda_Daughter_AntiPion_pz_DecayVertex	    = MomentumDaughterNegZ;
    Lambda_Daughter_AntiPion_pTPC		    = AntiPionTrack->GetTPCmomentum();
    Lambda_Daughter_AntiPion_Eta		    = AntiPionTrack->Eta();
    Lambda_Daughter_AntiPion_Phi		    = AntiPionTrack->Phi();
    Lambda_Daughter_AntiPion_TPC_Chi2		    = AntiPionTrack->GetTPCchi2();
    Lambda_Daughter_AntiPion_TPC_dEdx		    = AntiPionTrack->GetTPCsignal();
    Lambda_Daughter_AntiPion_TPC_dEdx_nSigma	    = (float)fPIDResponse->NumberOfSigmasTPC(AntiPionTrack,AliPID::kPion);
    Lambda_Daughter_AntiPion_TOF_Mass2		    = (float)AntiPion_TOF_m2;
    Lambda_Daughter_AntiPion_TOF_Mass2_nSigma	    = (float)AntiPion_TOF_m2_nSigma;
    Lambda_Daughter_AntiPion_ITS_dEdx		    = (float)AntiPion_ITS_dEdx;
    Lambda_Daughter_AntiPion_ITS_dEdx_nSigma	    = (float)AntiPion_ITS_dEdx_nSigma;
    Lambda_Daughter_AntiPion_DCAxy		    = AntiPion_DCAxy;
    Lambda_Daughter_AntiPion_DCAz		    = AntiPion_DCAz;
    Lambda_Daughter_AntiPion_TPC_nCrossedRows	    = AntiPionTrack->GetTPCCrossedRows();
    Lambda_Daughter_AntiPion_TPC_nSharedCluster	    = AntiPionTrack->GetTPCnclsS();
    Lambda_Daughter_AntiPion_TPC_nFindableCluster   = AntiPionTrack->GetTPCNclsF();
    Lambda_Daughter_AntiPion_TPC_nCluster	    = AntiPionTrack->GetTPCNcls();
    Lambda_Daughter_AntiPion_ITS_nCluster	    = (unsigned short)AntiPion_ITS_nCluster;
    Lambda_Daughter_AntiPion_ID			    = Label_Daughter2;


    fTempTree_Lambda->Fill();

    nLambdasSelected++;

  } // end of Loop loop




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ deuteron selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float     Deuteron_px;
  float     Deuteron_py;
  float     Deuteron_pz;
  float     Deuteron_px_Generated = 0.0;
  float     Deuteron_py_Generated = 0.0;
  float     Deuteron_pz_Generated = 0.0;
  float     Deuteron_pTPC;
  float     Deuteron_Eta;
  float     Deuteron_Phi;
  float     Deuteron_TPC_Chi2;
  float     Deuteron_TPC_dEdx;
  float     Deuteron_TPC_dEdx_nSigma;
  float     Deuteron_TOF_Mass2;
  float     Deuteron_TOF_Mass2_nSigma;
  float     Deuteron_ITS_dEdx;
  float     Deuteron_ITS_dEdx_nSigma;
  float     Deuteron_DCAxy;
  float     Deuteron_DCAz;
  unsigned short  Deuteron_TPC_nCrossedRows;
  unsigned short  Deuteron_TPC_nSharedCluster;
  unsigned short  Deuteron_TPC_nFindableCluster;
  unsigned short  Deuteron_TPC_nCluster;
  unsigned short  Deuteron_ITS_nCluster;
  int		  Deuteron_PDG;
  int		  Deuteron_MotherPDG;
  unsigned int	  Deuteron_ID;
  unsigned long   Deuteron_Event_Identifier;
  bool		  Deuteron_FilterBit = false;

  TTree *fTempTree_Deuteron = new TTree("fTempTree_Deuteron","fTempTree_Deuteron");
  fTempTree_Deuteron->Branch("Deuteron_px",&Deuteron_px,"Deuteron_px/F");
  fTempTree_Deuteron->Branch("Deuteron_py",&Deuteron_py,"Deuteron_py/F");
  fTempTree_Deuteron->Branch("Deuteron_pz",&Deuteron_pz,"Deuteron_pz/F");
  fTempTree_Deuteron->Branch("Deuteron_px_Generated",&Deuteron_px_Generated,"Deuteron_px_Generated/F");
  fTempTree_Deuteron->Branch("Deuteron_py_Generated",&Deuteron_py_Generated,"Deuteron_py_Generated/F");
  fTempTree_Deuteron->Branch("Deuteron_pz_Generated",&Deuteron_pz_Generated,"Deuteron_pz_Generated/F");
  fTempTree_Deuteron->Branch("Deuteron_pTPC",&Deuteron_pTPC,"Deuteron_pTPC/F");
  fTempTree_Deuteron->Branch("Deuteron_Eta",&Deuteron_Eta,"Deuteron_Eta/F");
  fTempTree_Deuteron->Branch("Deuteron_Phi",&Deuteron_Phi,"Deuteron_Phi/F");
  fTempTree_Deuteron->Branch("Deuteron_TPC_Chi2",&Deuteron_TPC_Chi2,"Deuteron_TPC_Chi2/F");
  fTempTree_Deuteron->Branch("Deuteron_TPC_dEdx",&Deuteron_TPC_dEdx,"Deuteron_TPC_dEdx/F");
  fTempTree_Deuteron->Branch("Deuteron_TPC_dEdx_nSigma",&Deuteron_TPC_dEdx_nSigma,"Deuteron_TPC_dEdx_nSigma/F");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Mass2",&Deuteron_TOF_Mass2,"Deuteron_TOF_Mass2/F");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Mass2_nSigma",&Deuteron_TOF_Mass2_nSigma,"Deuteron_TOF_Mass2_nSigma/F");
  fTempTree_Deuteron->Branch("Deuteron_ITS_dEdx",&Deuteron_ITS_dEdx,"Deuteron_ITS_dEdx/F");
  fTempTree_Deuteron->Branch("Deuteron_ITS_dEdx_nSigma",&Deuteron_ITS_dEdx_nSigma,"Deuteron_ITS_dEdx_nSigma/F");
  fTempTree_Deuteron->Branch("Deuteron_DCAxy",&Deuteron_DCAxy,"Deuteron_DCAxy/F");
  fTempTree_Deuteron->Branch("Deuteron_DCAz",&Deuteron_DCAz,"Deuteron_DCAz/F");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&Deuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&Deuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nFindableCluster",&Deuteron_TPC_nFindableCluster,"Deuteron_TPC_nFindableCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nCluster",&Deuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_ITS_nCluster",&Deuteron_ITS_nCluster,"Deuteron_ITS_nCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_PDG",&Deuteron_PDG,"Deuteron_PDG/I");
  fTempTree_Deuteron->Branch("Deuteron_MotherPDG",&Deuteron_MotherPDG,"Deuteron_MotherPDG/I");
  fTempTree_Deuteron->Branch("Deuteron_FilterBit",&Deuteron_FilterBit,"Deuteron_FilterBit/O");
  fTempTree_Deuteron->Branch("Deuteron_ID",&Deuteron_ID,"Deuteron_ID/i");
  fTempTree_Deuteron->Branch("Deuteron_Event_Identifier",&Deuteron_Event_Identifier,"Deuteron_Event_Identifier/l");

  unsigned short nDeuteronsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// deuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
  

    // apply deuteron cuts
    bool PassedDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,true,RunNumber);
    if(!PassedDeuteronCuts) continue;


    AliMCParticle *MCParticle = 0x0;
    int Label = TMath::Abs(Track->GetLabel());

    float Generated_px = 0.0;
    float Generated_py = 0.0;
    float Generated_pz = 0.0;
    int PDG = 0;
    int MotherPDG = 0;

    if(fIsMC == true){

      MCParticle = (AliMCParticle*) fMCEvent->GetTrack(Label);
      PDG = MCParticle->PdgCode();
      if(MCParticle->IsPhysicalPrimary() == true)         MotherPDG = 1;
      if(MCParticle->IsSecondaryFromMaterial() == true)   MotherPDG = 2;
      if(MCParticle->IsSecondaryFromWeakDecay() == true){
      
        int LabelMother = TMath::Abs(MCParticle->GetMother());
        AliMCParticle *MCParticleMother = (AliMCParticle*) fMCEvent->GetTrack(LabelMother);
        MotherPDG = MCParticleMother->PdgCode();

      } 

      Generated_px = MCParticle->Px();
      Generated_py = MCParticle->Py();
      Generated_pz = MCParticle->Pz();

    } // end of fIsMC == true






  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    float DCAxy = xv[0];
    float DCAz = xv[1];

    double TOF_m2	  = -999.0;
    double TOF_m2_nSigma  = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    bool TOFisOK = false;
    if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
    if(TOFisOK){

      TOF_m2	      = CalculateMassSquareTOF(*Track);
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,2,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,2,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();
      if(ITS_nCluster < 0) ITS_nCluster = 0;

    }


    Deuteron_px			    = Track->Px();
    Deuteron_py			    = Track->Py();
    Deuteron_pz			    = Track->Pz();
    Deuteron_px_Generated	    = Generated_px;
    Deuteron_py_Generated	    = Generated_py;
    Deuteron_pz_Generated	    = Generated_pz;
    Deuteron_pTPC		    = Track->GetTPCmomentum();
    Deuteron_Eta		    = Track->Eta();
    Deuteron_Phi		    = Track->Phi();
    Deuteron_TPC_Chi2		    = Track->GetTPCchi2();
    Deuteron_TPC_dEdx		    = Track->GetTPCsignal();
    Deuteron_TPC_dEdx_nSigma	    = (float)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    Deuteron_TOF_Mass2		    = (float)TOF_m2;
    Deuteron_TOF_Mass2_nSigma	    = (float)TOF_m2_nSigma;
    Deuteron_ITS_dEdx		    = (float)ITS_dEdx;
    Deuteron_ITS_dEdx_nSigma	    = (float)ITS_dEdx_nSigma;
    Deuteron_DCAxy		    = DCAxy;
    Deuteron_DCAz		    = DCAz;
    Deuteron_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    Deuteron_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    Deuteron_TPC_nFindableCluster   = Track->GetTPCNclsF();
    Deuteron_TPC_nCluster	    = Track->GetTPCNcls();
    Deuteron_ITS_nCluster	    = (unsigned short)ITS_nCluster;
    Deuteron_ID			    = track;
    Deuteron_Event_Identifier	    = EventID;
    Deuteron_PDG		    = PDG;
    Deuteron_MotherPDG		    = MotherPDG;
    Deuteron_FilterBit		    = Track->TestFilterBit(BIT(0));

    fTempTree_Deuteron->Fill();
    nDeuteronsSelected++;

  } // end of deuteron loop



  if((fSavePairsOnly == false) || ((fSavePairsOnly == true) && (nLambdasSelected > 0) && (nDeuteronsSelected > 0))){

    for(int Lambda = 0; Lambda < nLambdasSelected; Lambda++){

      TBranch *Branch_Lambda_px			  = fTempTree_Lambda->GetBranch("Lambda_px");
      TBranch *Branch_Lambda_py			  = fTempTree_Lambda->GetBranch("Lambda_py");
      TBranch *Branch_Lambda_pz			  = fTempTree_Lambda->GetBranch("Lambda_pz");
      TBranch *Branch_Lambda_px_Generated	  = fTempTree_Lambda->GetBranch("Lambda_px_Generated");
      TBranch *Branch_Lambda_py_Generated	  = fTempTree_Lambda->GetBranch("Lambda_py_Generated");
      TBranch *Branch_Lambda_pz_Generated	  = fTempTree_Lambda->GetBranch("Lambda_pz_Generated");
      TBranch *Branch_Lambda_Eta		  = fTempTree_Lambda->GetBranch("Lambda_Eta");
      TBranch *Branch_Lambda_Phi		  = fTempTree_Lambda->GetBranch("Lambda_Phi");
      TBranch *Branch_Lambda_TransverseRadius	  = fTempTree_Lambda->GetBranch("Lambda_TransverseRadius");
      TBranch *Branch_Lambda_CosinePointingAngle  = fTempTree_Lambda->GetBranch("Lambda_CosinePointingAngle");
      TBranch *Branch_Lambda_DCAv0ToPrimaryVertex = fTempTree_Lambda->GetBranch("Lambda_DCAv0ToPrimaryVertex");
      TBranch *Branch_Lambda_DCAv0Daughters	  = fTempTree_Lambda->GetBranch("Lambda_DCAv0Daughters");
      TBranch *Branch_Lambda_Alpha		  = fTempTree_Lambda->GetBranch("Lambda_Alpha");
      TBranch *Branch_Lambda_qT			  = fTempTree_Lambda->GetBranch("Lambda_qT");
      TBranch *Branch_Lambda_DecayLength	  = fTempTree_Lambda->GetBranch("Lambda_DecayLength");
      TBranch *Branch_Lambda_OpenAngle		  = fTempTree_Lambda->GetBranch("Lambda_OpenAngle");
      TBranch *Branch_Lambda_PDG_Daughter1	  = fTempTree_Lambda->GetBranch("Lambda_PDG_Daughter1");
      TBranch *Branch_Lambda_PDG_Daughter2	  = fTempTree_Lambda->GetBranch("Lambda_PDG_Daughter2");
      TBranch *Branch_Lambda_PDG_v01		  = fTempTree_Lambda->GetBranch("Lambda_PDG_v01");
      TBranch *Branch_Lambda_PDG_v02		  = fTempTree_Lambda->GetBranch("Lambda_PDG_v02");
      TBranch *Branch_Lambda_PDG_Mother1	  = fTempTree_Lambda->GetBranch("Lambda_PDG_Mother1");
      TBranch *Branch_Lambda_PDG_Mother2	  = fTempTree_Lambda->GetBranch("Lambda_PDG_Mother2");
      TBranch *Branch_Lambda_SameV0		  = fTempTree_Lambda->GetBranch("Lambda_SameV0");

      TBranch *Branch_Lambda_Daughter_Proton_px			  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_px");
      TBranch *Branch_Lambda_Daughter_Proton_py			  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_py");
      TBranch *Branch_Lambda_Daughter_Proton_pz			  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_pz");
      TBranch *Branch_Lambda_Daughter_Proton_px_Generated	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_px_Generated");
      TBranch *Branch_Lambda_Daughter_Proton_py_Generated	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_py_Generated");
      TBranch *Branch_Lambda_Daughter_Proton_pz_Generated	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_pz_Generated");
      TBranch *Branch_Lambda_Daughter_Proton_px_DecayVertex	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_px_DecayVertex");
      TBranch *Branch_Lambda_Daughter_Proton_py_DecayVertex	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_py_DecayVertex");
      TBranch *Branch_Lambda_Daughter_Proton_pz_DecayVertex	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_pz_DecayVertex");
      TBranch *Branch_Lambda_Daughter_Proton_pTPC		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_pTPC");
      TBranch *Branch_Lambda_Daughter_Proton_Eta		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_Eta");
      TBranch *Branch_Lambda_Daughter_Proton_Phi		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_Phi");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_Chi2		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_Chi2");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_dEdx		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_dEdx");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_dEdx_nSigma	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_dEdx_nSigma");
      TBranch *Branch_Lambda_Daughter_Proton_TOF_Mass2		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TOF_Mass2");
      TBranch *Branch_Lambda_Daughter_Proton_TOF_Mass2_nSigma	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TOF_Mass2_nSigma");
      TBranch *Branch_Lambda_Daughter_Proton_ITS_dEdx		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_ITS_dEdx");
      TBranch *Branch_Lambda_Daughter_Proton_ITS_dEdx_nSigma	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_ITS_dEdx_nSigma");
      TBranch *Branch_Lambda_Daughter_Proton_DCAxy		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_DCAxy");
      TBranch *Branch_Lambda_Daughter_Proton_DCAz		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_DCAz");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_nCrossedRows	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_nCrossedRows");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_nSharedCluster	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_nSharedCluster");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_nFindableCluster = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_nFindableCluster");
      TBranch *Branch_Lambda_Daughter_Proton_TPC_nCluster	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_TPC_nCluster");
      TBranch *Branch_Lambda_Daughter_Proton_ITS_nCluster	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_ITS_nCluster");
      TBranch *Branch_Lambda_Daughter_Proton_ID			  = fTempTree_Lambda->GetBranch("Lambda_Daughter_Proton_ID");

      TBranch *Branch_Lambda_Daughter_AntiPion_px		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_px");
      TBranch *Branch_Lambda_Daughter_AntiPion_py		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_py");
      TBranch *Branch_Lambda_Daughter_AntiPion_pz		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_pz");
      TBranch *Branch_Lambda_Daughter_AntiPion_px_Generated	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_px_Generated");
      TBranch *Branch_Lambda_Daughter_AntiPion_py_Generated	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_py_Generated");
      TBranch *Branch_Lambda_Daughter_AntiPion_pz_Generated	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_pz_Generated");
      TBranch *Branch_Lambda_Daughter_AntiPion_px_DecayVertex	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_px_DecayVertex");
      TBranch *Branch_Lambda_Daughter_AntiPion_py_DecayVertex	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_py_DecayVertex");
      TBranch *Branch_Lambda_Daughter_AntiPion_pz_DecayVertex	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_pz_DecayVertex");
      TBranch *Branch_Lambda_Daughter_AntiPion_pTPC		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_pTPC");
      TBranch *Branch_Lambda_Daughter_AntiPion_Eta		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_Eta");
      TBranch *Branch_Lambda_Daughter_AntiPion_Phi		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_Phi");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_Chi2		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_Chi2");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_dEdx		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_dEdx");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_dEdx_nSigma	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_dEdx_nSigma");
      TBranch *Branch_Lambda_Daughter_AntiPion_TOF_Mass2	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TOF_Mass2");
      TBranch *Branch_Lambda_Daughter_AntiPion_TOF_Mass2_nSigma	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TOF_Mass2_nSigma");
      TBranch *Branch_Lambda_Daughter_AntiPion_ITS_dEdx		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_ITS_dEdx");
      TBranch *Branch_Lambda_Daughter_AntiPion_ITS_dEdx_nSigma	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_ITS_dEdx_nSigma");
      TBranch *Branch_Lambda_Daughter_AntiPion_DCAxy		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_DCAxy");
      TBranch *Branch_Lambda_Daughter_AntiPion_DCAz		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_DCAz");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_nCrossedRows	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_nCrossedRows");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_nSharedCluster = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_nSharedCluster");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_nFindableCluster = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_nFindableCluster");
      TBranch *Branch_Lambda_Daughter_AntiPion_TPC_nCluster	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_TPC_nCluster");
      TBranch *Branch_Lambda_Daughter_AntiPion_ITS_nCluster	  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_ITS_nCluster");
      TBranch *Branch_Lambda_Daughter_AntiPion_ID		  = fTempTree_Lambda->GetBranch("Lambda_Daughter_AntiPion_ID");


      Branch_Lambda_px->SetAddress(&fLambda_px);
      Branch_Lambda_py->SetAddress(&fLambda_py);
      Branch_Lambda_pz->SetAddress(&fLambda_pz);
      Branch_Lambda_px_Generated->SetAddress(&fLambda_px_Generated);
      Branch_Lambda_py_Generated->SetAddress(&fLambda_py_Generated);
      Branch_Lambda_pz_Generated->SetAddress(&fLambda_pz_Generated);
      Branch_Lambda_Eta->SetAddress(&fLambda_Eta);
      Branch_Lambda_Phi->SetAddress(&fLambda_Phi);
      Branch_Lambda_TransverseRadius->SetAddress(&fLambda_TransverseRadius);
      Branch_Lambda_CosinePointingAngle->SetAddress(&fLambda_CosinePointingAngle);
      Branch_Lambda_DCAv0ToPrimaryVertex->SetAddress(&fLambda_DCAv0ToPrimaryVertex);
      Branch_Lambda_DCAv0Daughters->SetAddress(&fLambda_DCAv0Daughters);
      Branch_Lambda_Alpha->SetAddress(&fLambda_Alpha);
      Branch_Lambda_qT->SetAddress(&fLambda_qT);
      Branch_Lambda_DecayLength->SetAddress(&fLambda_DecayLength);
      Branch_Lambda_OpenAngle->SetAddress(&fLambda_OpenAngle);
      Branch_Lambda_PDG_Daughter1->SetAddress(&fLambda_PDG_Daughter1);
      Branch_Lambda_PDG_Daughter2->SetAddress(&fLambda_PDG_Daughter2);
      Branch_Lambda_PDG_v01->SetAddress(&fLambda_PDG_v01);
      Branch_Lambda_PDG_v02->SetAddress(&fLambda_PDG_v02);
      Branch_Lambda_PDG_Mother1->SetAddress(&fLambda_PDG_Mother1);
      Branch_Lambda_PDG_Mother2->SetAddress(&fLambda_PDG_Mother2);
      Branch_Lambda_SameV0->SetAddress(&fLambda_SameV0);

      Branch_Lambda_Daughter_Proton_px->SetAddress(&fLambda_Daughter_Proton_px);
      Branch_Lambda_Daughter_Proton_py->SetAddress(&fLambda_Daughter_Proton_py);
      Branch_Lambda_Daughter_Proton_pz->SetAddress(&fLambda_Daughter_Proton_pz);
      Branch_Lambda_Daughter_Proton_px_Generated->SetAddress(&fLambda_Daughter_Proton_px_Generated);
      Branch_Lambda_Daughter_Proton_py_Generated->SetAddress(&fLambda_Daughter_Proton_py_Generated);
      Branch_Lambda_Daughter_Proton_pz_Generated->SetAddress(&fLambda_Daughter_Proton_pz_Generated);
      Branch_Lambda_Daughter_Proton_px_DecayVertex->SetAddress(&fLambda_Daughter_Proton_px_DecayVertex);
      Branch_Lambda_Daughter_Proton_py_DecayVertex->SetAddress(&fLambda_Daughter_Proton_py_DecayVertex);
      Branch_Lambda_Daughter_Proton_pz_DecayVertex->SetAddress(&fLambda_Daughter_Proton_pz_DecayVertex);
      Branch_Lambda_Daughter_Proton_pTPC->SetAddress(&fLambda_Daughter_Proton_pTPC);
      Branch_Lambda_Daughter_Proton_Eta->SetAddress(&fLambda_Daughter_Proton_Eta);
      Branch_Lambda_Daughter_Proton_Phi->SetAddress(&fLambda_Daughter_Proton_Phi);
      Branch_Lambda_Daughter_Proton_TPC_Chi2->SetAddress(&fLambda_Daughter_Proton_TPC_Chi2);
      Branch_Lambda_Daughter_Proton_TPC_dEdx->SetAddress(&fLambda_Daughter_Proton_TPC_dEdx);
      Branch_Lambda_Daughter_Proton_TPC_dEdx_nSigma->SetAddress(&fLambda_Daughter_Proton_TPC_dEdx_nSigma);
      Branch_Lambda_Daughter_Proton_TOF_Mass2->SetAddress(&fLambda_Daughter_Proton_TOF_Mass2);
      Branch_Lambda_Daughter_Proton_TOF_Mass2_nSigma->SetAddress(&fLambda_Daughter_Proton_TOF_Mass2_nSigma);
      Branch_Lambda_Daughter_Proton_ITS_dEdx->SetAddress(&fLambda_Daughter_Proton_ITS_dEdx);
      Branch_Lambda_Daughter_Proton_ITS_dEdx_nSigma->SetAddress(&fLambda_Daughter_Proton_ITS_dEdx_nSigma);
      Branch_Lambda_Daughter_Proton_DCAxy->SetAddress(&fLambda_Daughter_Proton_DCAxy);
      Branch_Lambda_Daughter_Proton_DCAz->SetAddress(&fLambda_Daughter_Proton_DCAz);
      Branch_Lambda_Daughter_Proton_TPC_nCrossedRows->SetAddress(&fLambda_Daughter_Proton_TPC_nCrossedRows);
      Branch_Lambda_Daughter_Proton_TPC_nSharedCluster->SetAddress(&fLambda_Daughter_Proton_TPC_nSharedCluster);
      Branch_Lambda_Daughter_Proton_TPC_nFindableCluster->SetAddress(&fLambda_Daughter_Proton_TPC_nFindableCluster);
      Branch_Lambda_Daughter_Proton_TPC_nCluster->SetAddress(&fLambda_Daughter_Proton_TPC_nCluster);
      Branch_Lambda_Daughter_Proton_ITS_nCluster->SetAddress(&fLambda_Daughter_Proton_ITS_nCluster);
      Branch_Lambda_Daughter_Proton_ID->SetAddress(&fLambda_Daughter_Proton_ID);

      Branch_Lambda_Daughter_AntiPion_px->SetAddress(&fLambda_Daughter_AntiPion_px);
      Branch_Lambda_Daughter_AntiPion_py->SetAddress(&fLambda_Daughter_AntiPion_py);
      Branch_Lambda_Daughter_AntiPion_pz->SetAddress(&fLambda_Daughter_AntiPion_pz);
      Branch_Lambda_Daughter_AntiPion_px_Generated->SetAddress(&fLambda_Daughter_AntiPion_px_Generated);
      Branch_Lambda_Daughter_AntiPion_py_Generated->SetAddress(&fLambda_Daughter_AntiPion_py_Generated);
      Branch_Lambda_Daughter_AntiPion_pz_Generated->SetAddress(&fLambda_Daughter_AntiPion_pz_Generated);
      Branch_Lambda_Daughter_AntiPion_px_DecayVertex->SetAddress(&fLambda_Daughter_AntiPion_px_DecayVertex);
      Branch_Lambda_Daughter_AntiPion_py_DecayVertex->SetAddress(&fLambda_Daughter_AntiPion_py_DecayVertex);
      Branch_Lambda_Daughter_AntiPion_pz_DecayVertex->SetAddress(&fLambda_Daughter_AntiPion_pz_DecayVertex);
      Branch_Lambda_Daughter_AntiPion_pTPC->SetAddress(&fLambda_Daughter_AntiPion_pTPC);
      Branch_Lambda_Daughter_AntiPion_Eta->SetAddress(&fLambda_Daughter_AntiPion_Eta);
      Branch_Lambda_Daughter_AntiPion_Phi->SetAddress(&fLambda_Daughter_AntiPion_Phi);
      Branch_Lambda_Daughter_AntiPion_TPC_Chi2->SetAddress(&fLambda_Daughter_AntiPion_TPC_Chi2);
      Branch_Lambda_Daughter_AntiPion_TPC_dEdx->SetAddress(&fLambda_Daughter_AntiPion_TPC_dEdx);
      Branch_Lambda_Daughter_AntiPion_TPC_dEdx_nSigma->SetAddress(&fLambda_Daughter_AntiPion_TPC_dEdx_nSigma);
      Branch_Lambda_Daughter_AntiPion_TOF_Mass2->SetAddress(&fLambda_Daughter_AntiPion_TOF_Mass2);
      Branch_Lambda_Daughter_AntiPion_TOF_Mass2_nSigma->SetAddress(&fLambda_Daughter_AntiPion_TOF_Mass2_nSigma);
      Branch_Lambda_Daughter_AntiPion_ITS_dEdx->SetAddress(&fLambda_Daughter_AntiPion_ITS_dEdx);
      Branch_Lambda_Daughter_AntiPion_ITS_dEdx_nSigma->SetAddress(&fLambda_Daughter_AntiPion_ITS_dEdx_nSigma);
      Branch_Lambda_Daughter_AntiPion_DCAxy->SetAddress(&fLambda_Daughter_AntiPion_DCAxy);
      Branch_Lambda_Daughter_AntiPion_DCAz->SetAddress(&fLambda_Daughter_AntiPion_DCAz);
      Branch_Lambda_Daughter_AntiPion_TPC_nCrossedRows->SetAddress(&fLambda_Daughter_AntiPion_TPC_nCrossedRows);
      Branch_Lambda_Daughter_AntiPion_TPC_nSharedCluster->SetAddress(&fLambda_Daughter_AntiPion_TPC_nSharedCluster);
      Branch_Lambda_Daughter_AntiPion_TPC_nFindableCluster->SetAddress(&fLambda_Daughter_AntiPion_TPC_nFindableCluster);
      Branch_Lambda_Daughter_AntiPion_TPC_nCluster->SetAddress(&fLambda_Daughter_AntiPion_TPC_nCluster);
      Branch_Lambda_Daughter_AntiPion_ITS_nCluster->SetAddress(&fLambda_Daughter_AntiPion_ITS_nCluster);
      Branch_Lambda_Daughter_AntiPion_ID->SetAddress(&fLambda_Daughter_AntiPion_ID);



      Branch_Lambda_px->SetAutoDelete(true);
      Branch_Lambda_py->SetAutoDelete(true);
      Branch_Lambda_pz->SetAutoDelete(true);
      Branch_Lambda_px_Generated->SetAutoDelete(true);
      Branch_Lambda_py_Generated->SetAutoDelete(true);
      Branch_Lambda_pz_Generated->SetAutoDelete(true);
      Branch_Lambda_Eta->SetAutoDelete(true);
      Branch_Lambda_Phi->SetAutoDelete(true);
      Branch_Lambda_TransverseRadius->SetAutoDelete(true);
      Branch_Lambda_CosinePointingAngle->SetAutoDelete(true);
      Branch_Lambda_DCAv0ToPrimaryVertex->SetAutoDelete(true);
      Branch_Lambda_DCAv0Daughters->SetAutoDelete(true);
      Branch_Lambda_Alpha->SetAutoDelete(true);
      Branch_Lambda_qT->SetAutoDelete(true);
      Branch_Lambda_DecayLength->SetAutoDelete(true);
      Branch_Lambda_OpenAngle->SetAutoDelete(true);
      Branch_Lambda_PDG_Daughter1->SetAutoDelete(true);
      Branch_Lambda_PDG_Daughter2->SetAutoDelete(true);
      Branch_Lambda_PDG_v01->SetAutoDelete(true);
      Branch_Lambda_PDG_v02->SetAutoDelete(true);
      Branch_Lambda_PDG_Mother1->SetAutoDelete(true);
      Branch_Lambda_PDG_Mother2->SetAutoDelete(true);
      Branch_Lambda_SameV0->SetAutoDelete(true);

      Branch_Lambda_Daughter_Proton_px->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_py->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_pz->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_px_Generated->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_py_Generated->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_pz_Generated->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_px_DecayVertex->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_py_DecayVertex->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_pz_DecayVertex->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_pTPC->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_Eta->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_Phi->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_Chi2->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_dEdx->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TOF_Mass2->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_ITS_dEdx->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_DCAxy->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_DCAz->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_nFindableCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_TPC_nCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_ITS_nCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_Proton_ID->SetAutoDelete(true);

      Branch_Lambda_Daughter_AntiPion_px->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_py->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_pz->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_px_Generated->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_py_Generated->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_pz_Generated->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_px_DecayVertex->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_py_DecayVertex->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_pz_DecayVertex->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_pTPC->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_Eta->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_Phi->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_Chi2->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_dEdx->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TOF_Mass2->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_ITS_dEdx->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_DCAxy->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_DCAz->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_nFindableCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_TPC_nCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_ITS_nCluster->SetAutoDelete(true);
      Branch_Lambda_Daughter_AntiPion_ID->SetAutoDelete(true);

      fTempTree_Lambda->GetEntry(Lambda);

      fLambda_Event_Multiplicity    = Multiplicity;
      fLambda_Event_Centrality	    = Centrality;
      fLambda_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fLambda_Event_BField	    = BField;
      fLambda_Event_Identifier	    = EventID;
      fLambda_Event_ParticleNumber  = Lambda;
      fLambda_Event_nParticles	    = nLambdasSelected;

      fSaveTree_Lambda->Fill();

    } // end of loop (copy Lambdas)


    for(int Deuteron = 0; Deuteron < nDeuteronsSelected; Deuteron++){

      TBranch *Branch_Deuteron_px		    = fTempTree_Deuteron->GetBranch("Deuteron_px");
      TBranch *Branch_Deuteron_py		    = fTempTree_Deuteron->GetBranch("Deuteron_py");
      TBranch *Branch_Deuteron_pz		    = fTempTree_Deuteron->GetBranch("Deuteron_pz");
      TBranch *Branch_Deuteron_px_Generated	    = fTempTree_Deuteron->GetBranch("Deuteron_px_Generated");
      TBranch *Branch_Deuteron_py_Generated	    = fTempTree_Deuteron->GetBranch("Deuteron_py_Generated");
      TBranch *Branch_Deuteron_pz_Generated	    = fTempTree_Deuteron->GetBranch("Deuteron_pz_Generated");
      TBranch *Branch_Deuteron_pTPC		    = fTempTree_Deuteron->GetBranch("Deuteron_pTPC");
      TBranch *Branch_Deuteron_Eta		    = fTempTree_Deuteron->GetBranch("Deuteron_Eta");
      TBranch *Branch_Deuteron_Phi		    = fTempTree_Deuteron->GetBranch("Deuteron_Phi");
      TBranch *Branch_Deuteron_TPC_Chi2		    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_Chi2");
      TBranch *Branch_Deuteron_TPC_dEdx		    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_dEdx");
      TBranch *Branch_Deuteron_TPC_dEdx_nSigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_dEdx_nSigma");
      TBranch *Branch_Deuteron_TOF_Mass2	    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Mass2");
      TBranch *Branch_Deuteron_TOF_Mass2_nSigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Mass2_nSigma");
      TBranch *Branch_Deuteron_ITS_dEdx		    = fTempTree_Deuteron->GetBranch("Deuteron_ITS_dEdx");
      TBranch *Branch_Deuteron_ITS_dEdx_nSigma	    = fTempTree_Deuteron->GetBranch("Deuteron_ITS_dEdx_nSigma");
      TBranch *Branch_Deuteron_DCAxy		    = fTempTree_Deuteron->GetBranch("Deuteron_DCAxy");
      TBranch *Branch_Deuteron_DCAz		    = fTempTree_Deuteron->GetBranch("Deuteron_DCAz");
      TBranch *Branch_Deuteron_TPC_nCrossedRows	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nCrossedRows");
      TBranch *Branch_Deuteron_TPC_nSharedCluster   = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nSharedCluster");
      TBranch *Branch_Deuteron_TPC_nFindableCluster = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nFindableCluster");
      TBranch *Branch_Deuteron_TPC_nCluster	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nCluster");
      TBranch *Branch_Deuteron_ITS_nCluster	    = fTempTree_Deuteron->GetBranch("Deuteron_ITS_nCluster");
      TBranch *Branch_Deuteron_PDG		    = fTempTree_Deuteron->GetBranch("Deuteron_PDG");
      TBranch *Branch_Deuteron_MotherPDG	    = fTempTree_Deuteron->GetBranch("Deuteron_MotherPDG");
      TBranch *Branch_Deuteron_FilterBit	    = fTempTree_Deuteron->GetBranch("Deuteron_FilterBit");
      TBranch *Branch_Deuteron_ID		    = fTempTree_Deuteron->GetBranch("Deuteron_ID");
      TBranch *Branch_Deuteron_Event_Identifier	    = fTempTree_Deuteron->GetBranch("Deuteron_Event_Identifier");


      Branch_Deuteron_px->SetAddress(&fDeuteron_px);
      Branch_Deuteron_py->SetAddress(&fDeuteron_py);
      Branch_Deuteron_pz->SetAddress(&fDeuteron_pz);
      Branch_Deuteron_px_Generated->SetAddress(&fDeuteron_px_Generated);
      Branch_Deuteron_py_Generated->SetAddress(&fDeuteron_py_Generated);
      Branch_Deuteron_pz_Generated->SetAddress(&fDeuteron_pz_Generated);
      Branch_Deuteron_pTPC->SetAddress(&fDeuteron_pTPC);
      Branch_Deuteron_Eta->SetAddress(&fDeuteron_Eta);
      Branch_Deuteron_Phi->SetAddress(&fDeuteron_Phi);
      Branch_Deuteron_TPC_Chi2->SetAddress(&fDeuteron_TPC_Chi2);
      Branch_Deuteron_TPC_dEdx->SetAddress(&fDeuteron_TPC_dEdx);
      Branch_Deuteron_TPC_dEdx_nSigma->SetAddress(&fDeuteron_TPC_dEdx_nSigma);
      Branch_Deuteron_TOF_Mass2->SetAddress(&fDeuteron_TOF_Mass2);
      Branch_Deuteron_TOF_Mass2_nSigma->SetAddress(&fDeuteron_TOF_Mass2_nSigma);
      Branch_Deuteron_ITS_dEdx->SetAddress(&fDeuteron_ITS_dEdx);
      Branch_Deuteron_ITS_dEdx_nSigma->SetAddress(&fDeuteron_ITS_dEdx_nSigma);
      Branch_Deuteron_DCAxy->SetAddress(&fDeuteron_DCAxy);
      Branch_Deuteron_DCAz->SetAddress(&fDeuteron_DCAz);
      Branch_Deuteron_TPC_nCrossedRows->SetAddress(&fDeuteron_TPC_nCrossedRows);
      Branch_Deuteron_TPC_nSharedCluster->SetAddress(&fDeuteron_TPC_nSharedCluster);
      Branch_Deuteron_TPC_nFindableCluster->SetAddress(&fDeuteron_TPC_nFindableCluster);
      Branch_Deuteron_TPC_nCluster->SetAddress(&fDeuteron_TPC_nCluster);
      Branch_Deuteron_ITS_nCluster->SetAddress(&fDeuteron_ITS_nCluster);
      Branch_Deuteron_PDG->SetAddress(&fDeuteron_PDG);
      Branch_Deuteron_MotherPDG->SetAddress(&fDeuteron_MotherPDG);
      Branch_Deuteron_FilterBit->SetAddress(&fDeuteron_FilterBit);
      Branch_Deuteron_ID->SetAddress(&fDeuteron_ID);
      Branch_Deuteron_Event_Identifier->SetAddress(&fDeuteron_Event_Identifier);


      Branch_Deuteron_px->SetAutoDelete(true);
      Branch_Deuteron_py->SetAutoDelete(true);
      Branch_Deuteron_pz->SetAutoDelete(true);
      Branch_Deuteron_pTPC->SetAutoDelete(true);
      Branch_Deuteron_Eta->SetAutoDelete(true);
      Branch_Deuteron_Phi->SetAutoDelete(true);
      Branch_Deuteron_TPC_Chi2->SetAutoDelete(true);
      Branch_Deuteron_TPC_dEdx->SetAutoDelete(true);
      Branch_Deuteron_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_Deuteron_TOF_Mass2->SetAutoDelete(true);
      Branch_Deuteron_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_Deuteron_ITS_dEdx->SetAutoDelete(true);
      Branch_Deuteron_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_Deuteron_DCAxy->SetAutoDelete(true);
      Branch_Deuteron_DCAz->SetAutoDelete(true);
      Branch_Deuteron_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_Deuteron_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_Deuteron_TPC_nFindableCluster->SetAutoDelete(true);
      Branch_Deuteron_TPC_nCluster->SetAutoDelete(true);
      Branch_Deuteron_ITS_nCluster->SetAutoDelete(true);
      Branch_Deuteron_PDG->SetAutoDelete(true);
      Branch_Deuteron_MotherPDG->SetAutoDelete(true);
      Branch_Deuteron_FilterBit->SetAutoDelete(true);
      Branch_Deuteron_ID->SetAutoDelete(true);
      Branch_Deuteron_Event_Identifier->SetAutoDelete(true);

      fTempTree_Deuteron->GetEntry(Deuteron);


      fDeuteron_Event_Multiplicity    = Multiplicity;
      fDeuteron_Event_Centrality      = Centrality;
      fDeuteron_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fDeuteron_Event_BField	      = BField;
      fDeuteron_Event_ParticleNumber  = Deuteron;
      fDeuteron_Event_nParticles      = nDeuteronsSelected;

      fSaveTree_Deuteron->Fill();


    } // end of loop (copy deuterons)


  } // end of same-event

  fTempTree_Lambda->Delete();
  fTempTree_Deuteron->Delete();








  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ ++ AntiLambda selection loop ++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    float AntiLambda_px;
    float AntiLambda_py;
    float AntiLambda_pz;
    float AntiLambda_px_Generated = 0.0;
    float AntiLambda_py_Generated = 0.0;
    float AntiLambda_pz_Generated = 0.0;
    float AntiLambda_Eta;
    float AntiLambda_Phi;
    float AntiLambda_TransverseRadius;
    float AntiLambda_CosinePointingAngle;
    float AntiLambda_DCAv0ToPrimaryVertex;
    float AntiLambda_DCAv0Daughters;
    float AntiLambda_Alpha;
    float AntiLambda_qT;
    float AntiLambda_DecayLength;
    float AntiLambda_OpenAngle;
    int	  AntiLambda_PDG_Daughter1;
    int	  AntiLambda_PDG_Daughter2;
    int	  AntiLambda_PDG_v01;
    int	  AntiLambda_PDG_v02;
    int	  AntiLambda_PDG_Mother1;
    int	  AntiLambda_PDG_Mother2;
    bool  AntiLambda_SameV0;

  float     AntiLambda_Daughter_AntiProton_px;
  float     AntiLambda_Daughter_AntiProton_py;
  float     AntiLambda_Daughter_AntiProton_pz;
  float     AntiLambda_Daughter_AntiProton_px_Generated;
  float     AntiLambda_Daughter_AntiProton_py_Generated;
  float     AntiLambda_Daughter_AntiProton_pz_Generated;
  float     AntiLambda_Daughter_AntiProton_px_DecayVertex;
  float     AntiLambda_Daughter_AntiProton_py_DecayVertex;
  float     AntiLambda_Daughter_AntiProton_pz_DecayVertex;
  float     AntiLambda_Daughter_AntiProton_pTPC;
  float     AntiLambda_Daughter_AntiProton_Eta;
  float     AntiLambda_Daughter_AntiProton_Phi;
  float     AntiLambda_Daughter_AntiProton_TPC_Chi2;
  float     AntiLambda_Daughter_AntiProton_TPC_dEdx;
  float     AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma;
  float     AntiLambda_Daughter_AntiProton_TOF_Mass2;
  float     AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma;
  float     AntiLambda_Daughter_AntiProton_ITS_dEdx;
  float     AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma;
  float     AntiLambda_Daughter_AntiProton_DCAxy;
  float     AntiLambda_Daughter_AntiProton_DCAz;
  unsigned short    AntiLambda_Daughter_AntiProton_TPC_nCrossedRows;
  unsigned short    AntiLambda_Daughter_AntiProton_TPC_nSharedCluster;
  unsigned short    AntiLambda_Daughter_AntiProton_TPC_nFindableCluster;
  unsigned short    AntiLambda_Daughter_AntiProton_TPC_nCluster;
  unsigned short    AntiLambda_Daughter_AntiProton_ITS_nCluster;
  unsigned int      AntiLambda_Daughter_AntiProton_ID = 1;

  float     AntiLambda_Daughter_Pion_px;
  float     AntiLambda_Daughter_Pion_py;
  float     AntiLambda_Daughter_Pion_pz;
  float     AntiLambda_Daughter_Pion_px_Generated;
  float     AntiLambda_Daughter_Pion_py_Generated;
  float     AntiLambda_Daughter_Pion_pz_Generated;
  float     AntiLambda_Daughter_Pion_px_DecayVertex;
  float     AntiLambda_Daughter_Pion_py_DecayVertex;
  float     AntiLambda_Daughter_Pion_pz_DecayVertex;
  float     AntiLambda_Daughter_Pion_pTPC;
  float     AntiLambda_Daughter_Pion_Eta;
  float     AntiLambda_Daughter_Pion_Phi;
  float     AntiLambda_Daughter_Pion_TPC_Chi2;
  float     AntiLambda_Daughter_Pion_TPC_dEdx;
  float     AntiLambda_Daughter_Pion_TPC_dEdx_nSigma;
  float     AntiLambda_Daughter_Pion_TOF_Mass2;
  float     AntiLambda_Daughter_Pion_TOF_Mass2_nSigma;
  float     AntiLambda_Daughter_Pion_ITS_dEdx;
  float     AntiLambda_Daughter_Pion_ITS_dEdx_nSigma;
  float     AntiLambda_Daughter_Pion_DCAxy;
  float     AntiLambda_Daughter_Pion_DCAz;
  unsigned short    AntiLambda_Daughter_Pion_TPC_nCrossedRows;
  unsigned short    AntiLambda_Daughter_Pion_TPC_nSharedCluster;
  unsigned short    AntiLambda_Daughter_Pion_TPC_nFindableCluster;
  unsigned short    AntiLambda_Daughter_Pion_TPC_nCluster;
  unsigned short    AntiLambda_Daughter_Pion_ITS_nCluster;
  unsigned int      AntiLambda_Daughter_Pion_ID = 1;


  TTree *fTempTree_AntiLambda = new TTree("fTempTree_AntiLambda","fTempTree_AntiLambda");
  fTempTree_AntiLambda->Branch("AntiLambda_px",&AntiLambda_px,"AntiLambda_px/F");
  fTempTree_AntiLambda->Branch("AntiLambda_py",&AntiLambda_py,"AntiLambda_py/F");
  fTempTree_AntiLambda->Branch("AntiLambda_pz",&AntiLambda_pz,"AntiLambda_pz/F");
  fTempTree_AntiLambda->Branch("AntiLambda_px_Generated",&AntiLambda_px_Generated,"AntiLambda_px_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_py_Generated",&AntiLambda_py_Generated,"AntiLambda_py_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_pz_Generated",&AntiLambda_pz_Generated,"AntiLambda_pz_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Eta",&AntiLambda_Eta,"AntiLambda_Eta/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Phi",&AntiLambda_Phi,"AntiLambda_Phi/F");
  fTempTree_AntiLambda->Branch("AntiLambda_TransverseRadius",&AntiLambda_TransverseRadius,"AntiLambda_TransverseRadius/F");
  fTempTree_AntiLambda->Branch("AntiLambda_CosinePointingAngle",&AntiLambda_CosinePointingAngle,"AntiLambda_CosinePointingAngle/F");
  fTempTree_AntiLambda->Branch("AntiLambda_DCAv0ToPrimaryVertex",&AntiLambda_DCAv0ToPrimaryVertex,"AntiLambda_DCAv0ToPrimaryVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_DCAv0Daughters",&AntiLambda_DCAv0Daughters,"AntiLambda_DCAv0Daughters/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Alpha",&AntiLambda_Alpha,"AntiLambda_Alpha/F");
  fTempTree_AntiLambda->Branch("AntiLambda_qT",&AntiLambda_qT,"AntiLambda_qT/F");
  fTempTree_AntiLambda->Branch("AntiLambda_DecayLength",&AntiLambda_DecayLength,"AntiLambda_DecayLength/F");
  fTempTree_AntiLambda->Branch("AntiLambda_OpenAngle",&AntiLambda_OpenAngle,"AntiLambda_OpenAngle/F");
  fTempTree_AntiLambda->Branch("AntiLambda_PDG_Daughter1",&AntiLambda_PDG_Daughter1,"AntiLambda_PDG_Daughter1/I");
  fTempTree_AntiLambda->Branch("AntiLambda_PDG_Daughter2",&AntiLambda_PDG_Daughter2,"AntiLambda_PDG_Daughter2/I");
  fTempTree_AntiLambda->Branch("AntiLambda_PDG_v01",&AntiLambda_PDG_v01,"AntiLambda_PDG_v01/I");
  fTempTree_AntiLambda->Branch("AntiLambda_PDG_v02",&AntiLambda_PDG_v02,"AntiLambda_PDG_v02/I");
  fTempTree_AntiLambda->Branch("AntiLambda_PDG_Mother1",&AntiLambda_PDG_Mother1,"AntiLambda_PDG_Mother1/I");
  fTempTree_AntiLambda->Branch("AntiLambda_PDG_Mother2",&AntiLambda_PDG_Mother2,"AntiLambda_PDG_Mother2/I");
  fTempTree_AntiLambda->Branch("AntiLambda_SameV0",&AntiLambda_SameV0,"AntiLambda_SameV0/O");

  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px",&AntiLambda_Daughter_AntiProton_px,"AntiLambda_Daughter_AntiProton_px/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py",&AntiLambda_Daughter_AntiProton_py,"AntiLambda_Daughter_AntiProton_py/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz",&AntiLambda_Daughter_AntiProton_pz,"AntiLambda_Daughter_AntiProton_pz/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px_Generated",&AntiLambda_Daughter_AntiProton_px_Generated,"AntiLambda_Daughter_AntiProton_px_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py_Generated",&AntiLambda_Daughter_AntiProton_py_Generated,"AntiLambda_Daughter_AntiProton_py_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz_Generated",&AntiLambda_Daughter_AntiProton_pz_Generated,"AntiLambda_Daughter_AntiProton_pz_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px_DecayVertex",&AntiLambda_Daughter_AntiProton_px_DecayVertex,"AntiLambda_Daughter_AntiProton_px_DecayVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py_DecayVertex",&AntiLambda_Daughter_AntiProton_py_DecayVertex,"AntiLambda_Daughter_AntiProton_py_DecayVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz_DecayVertex",&AntiLambda_Daughter_AntiProton_pz_DecayVertex,"AntiLambda_Daughter_AntiProton_pz_DecayVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pTPC",&AntiLambda_Daughter_AntiProton_pTPC,"AntiLambda_Daughter_AntiProton_pTPC/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_Eta",&AntiLambda_Daughter_AntiProton_Eta,"AntiLambda_Daughter_AntiProton_Eta/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_Phi",&AntiLambda_Daughter_AntiProton_Phi,"AntiLambda_Daughter_AntiProton_Phi/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_Chi2",&AntiLambda_Daughter_AntiProton_TPC_Chi2,"AntiLambda_Daughter_AntiProton_TPC_Chi2/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_dEdx",&AntiLambda_Daughter_AntiProton_TPC_dEdx,"AntiLambda_Daughter_AntiProton_TPC_dEdx/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma",&AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma,"AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TOF_Mass2",&AntiLambda_Daughter_AntiProton_TOF_Mass2,"AntiLambda_Daughter_AntiProton_TOF_Mass2/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma",&AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma,"AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_dEdx",&AntiLambda_Daughter_AntiProton_ITS_dEdx,"AntiLambda_Daughter_AntiProton_ITS_dEdx/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma",&AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma,"AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_DCAxy",&AntiLambda_Daughter_AntiProton_DCAxy,"AntiLambda_Daughter_AntiProton_DCAxy/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_DCAz",&AntiLambda_Daughter_AntiProton_DCAz,"AntiLambda_Daughter_AntiProton_DCAz/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nCrossedRows",&AntiLambda_Daughter_AntiProton_TPC_nCrossedRows,"AntiLambda_Daughter_AntiProton_TPC_nCrossedRows/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nSharedCluster",&AntiLambda_Daughter_AntiProton_TPC_nSharedCluster,"AntiLambda_Daughter_AntiProton_TPC_nSharedCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nFindableCluster",&AntiLambda_Daughter_AntiProton_TPC_nFindableCluster,"AntiLambda_Daughter_AntiProton_TPC_nFindableCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nCluster",&AntiLambda_Daughter_AntiProton_TPC_nCluster,"AntiLambda_Daughter_AntiProton_TPC_nCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_nCluster",&AntiLambda_Daughter_AntiProton_ITS_nCluster,"AntiLambda_Daughter_AntiProton_ITS_nCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ID",&AntiLambda_Daughter_AntiProton_ID,"AntiLambda_Daughter_AntiProton_ID/i");

  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px",&AntiLambda_Daughter_Pion_px,"AntiLambda_Daughter_Pion_px/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py",&AntiLambda_Daughter_Pion_py,"AntiLambda_Daughter_Pion_py/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz",&AntiLambda_Daughter_Pion_pz,"AntiLambda_Daughter_Pion_pz/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px_Generated",&AntiLambda_Daughter_Pion_px_Generated,"AntiLambda_Daughter_Pion_px_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py_Generated",&AntiLambda_Daughter_Pion_py_Generated,"AntiLambda_Daughter_Pion_py_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz_Generated",&AntiLambda_Daughter_Pion_pz_Generated,"AntiLambda_Daughter_Pion_pz_Generated/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px_DecayVertex",&AntiLambda_Daughter_Pion_px_DecayVertex,"AntiLambda_Daughter_Pion_px_DecayVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py_DecayVertex",&AntiLambda_Daughter_Pion_py_DecayVertex,"AntiLambda_Daughter_Pion_py_DecayVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz_DecayVertex",&AntiLambda_Daughter_Pion_pz_DecayVertex,"AntiLambda_Daughter_Pion_pz_DecayVertex/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pTPC",&AntiLambda_Daughter_Pion_pTPC,"AntiLambda_Daughter_Pion_pTPC/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_Eta",&AntiLambda_Daughter_Pion_Eta,"AntiLambda_Daughter_Pion_Eta/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_Phi",&AntiLambda_Daughter_Pion_Phi,"AntiLambda_Daughter_Pion_Phi/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_Chi2",&AntiLambda_Daughter_Pion_TPC_Chi2,"AntiLambda_Daughter_Pion_TPC_Chi2/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_dEdx",&AntiLambda_Daughter_Pion_TPC_dEdx,"AntiLambda_Daughter_Pion_TPC_dEdx/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_dEdx_nSigma",&AntiLambda_Daughter_Pion_TPC_dEdx_nSigma,"AntiLambda_Daughter_Pion_TPC_dEdx_nSigma/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TOF_Mass2",&AntiLambda_Daughter_Pion_TOF_Mass2,"AntiLambda_Daughter_Pion_TOF_Mass2/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TOF_Mass2_nSigma",&AntiLambda_Daughter_Pion_TOF_Mass2_nSigma,"AntiLambda_Daughter_Pion_TOF_Mass2_nSigma/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_dEdx",&AntiLambda_Daughter_Pion_ITS_dEdx,"AntiLambda_Daughter_Pion_ITS_dEdx/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_dEdx_nSigma",&AntiLambda_Daughter_Pion_ITS_dEdx_nSigma,"AntiLambda_Daughter_Pion_ITS_dEdx_nSigma/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_DCAxy",&AntiLambda_Daughter_Pion_DCAxy,"AntiLambda_Daughter_Pion_DCAxy/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_DCAz",&AntiLambda_Daughter_Pion_DCAz,"AntiLambda_Daughter_Pion_DCAz/F");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nCrossedRows",&AntiLambda_Daughter_Pion_TPC_nCrossedRows,"AntiLambda_Daughter_Pion_TPC_nCrossedRows/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nSharedCluster",&AntiLambda_Daughter_Pion_TPC_nSharedCluster,"AntiLambda_Daughter_Pion_TPC_nSharedCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nFindableCluster",&AntiLambda_Daughter_Pion_TPC_nFindableCluster,"AntiLambda_Daughter_Pion_TPC_nFindableCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nCluster",&AntiLambda_Daughter_Pion_TPC_nCluster,"AntiLambda_Daughter_Pion_TPC_nCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_nCluster",&AntiLambda_Daughter_Pion_ITS_nCluster,"AntiLambda_Daughter_Pion_ITS_nCluster/s");
  fTempTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ID",&AntiLambda_Daughter_Pion_ID,"AntiLambda_Daughter_Pion_ID/i");





  unsigned short nAntiLambdasSelected = 0;

  for(int V0 = 0; V0 < nV0s; V0++)  // AntiLambda loop
  { 

    AliAODv0 *v0 = dynamic_cast<AliAODv0*>(fAODEvent->GetV0(V0));
    if(!v0)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODv0 found");
    if(!v0) continue;

    AliAODTrack *AntiProtonTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    if(!AntiProtonTrack) continue;

    AliAODTrack *PionTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    if(!PionTrack) continue;
 

    bool PassedAntiProtonCuts = CheckProtonCuts(*AntiProtonTrack,*fPIDResponse,false,RunNumber);
    if(!PassedAntiProtonCuts) continue;

    double PionDCAxy_min = FitPionDCAxy_min->Eval(PionTrack->Pt());
    double PionDCAz_min = FitPionDCAz_min->Eval(PionTrack->Pt());
    bool PassedPionCuts = CheckPionCuts(*PionTrack,*fPIDResponse,true,RunNumber,PionDCAxy_min,PionDCAz_min);
    if(!PassedPionCuts) continue;


    AliAODVertex *DecayVertex = v0->GetSecondaryVtx();

    AliExternalTrackParam *AntiProtonTrackParam = new AliExternalTrackParam();
    AntiProtonTrackParam->CopyFromVTrack(AntiProtonTrack);
    double AntiProton_DCA_new[2];
    double AntiProton_CovarianceMatrix_new[3];
    AntiProtonTrackParam->PropagateToDCA(DecayVertex,BField,10.0,AntiProton_DCA_new,AntiProton_CovarianceMatrix_new);
    double AntiProtonMomentumPropagated[3];
    AntiProtonTrackParam->GetPxPyPz(AntiProtonMomentumPropagated);

    AliExternalTrackParam *PionTrackParam = new AliExternalTrackParam();
    PionTrackParam->CopyFromVTrack(PionTrack);
    double Pion_DCA_new[2];
    double Pion_CovarianceMatrix_new[3];
    PionTrackParam->PropagateToDCA(DecayVertex,BField,10.0,Pion_DCA_new,Pion_CovarianceMatrix_new);
    double PionMomentumPropagated[3];
    PionTrackParam->GetPxPyPz(PionMomentumPropagated);


    bool PrintDecayMomentaOnScreen = false;
    if(PrintDecayMomentaOnScreen == true){

      std::cout << "\nAntiProtonTrack->Px() = " << AntiProtonTrack->Px() << "\t AntiProtonTrack->Py() = " << AntiProtonTrack->Py() << "\t AntiProtonTrack->Pz() = " << AntiProtonTrack->Pz() << std::endl;
      std::cout << "v0->MomNegX() = " << v0->MomNegX() << "\t v0->MomNegY() = " << v0->MomNegY() << "\tv0->MomNegZ() = " << v0->MomNegZ() << std::endl;
      std::cout << "px_prop = " << AntiProtonMomentumPropagated[0] << "\t\t py_prop = " << AntiProtonMomentumPropagated[1] << "\t\t pz_prop = " << AntiProtonMomentumPropagated[2] << std::endl;
    std::cout << "PionTrack->Px() = " << PionTrack->Px() << "\t PionTrack->Py() = " << PionTrack->Py() << "\tPionTrack->Pz() = " << PionTrack->Pz() << std::endl;
    std::cout << "v0->MomPosX() = " << v0->MomPosX() << "\t v0->MomPosY() = " << v0->MomPosY() << "\tv0->MomPosZ() = " << v0->MomPosZ() << std::endl;
    std::cout << "px_prop = " << PionMomentumPropagated[0] << "\t py_prop = " << PionMomentumPropagated[1] << "\t pz_prop = " << PionMomentumPropagated[2] << std::endl;

    } // end of PrintDecayMomentaOnScreen


    double MomentumDaughterNegX = AntiProtonMomentumPropagated[0];
    double MomentumDaughterNegY = AntiProtonMomentumPropagated[1];
    double MomentumDaughterNegZ = AntiProtonMomentumPropagated[2];
    double MomentumDaughterPosX = PionMomentumPropagated[0];
    double MomentumDaughterPosY = PionMomentumPropagated[1];
    double MomentumDaughterPosZ = PionMomentumPropagated[2];
    double MomentumV0X = v0->MomV0X();
    double MomentumV0Y = v0->MomV0Y();
    double MomentumV0Z = v0->MomV0Z();

    // apply AntiLambda cuts
    double MassInvariant = CalculateInvariantMassLambda(PionMomentumPropagated,AntiProtonMomentumPropagated,false);
    double WrongMassInvariant = CalculateInvariantMassLambda(AntiProtonMomentumPropagated,PionMomentumPropagated,false);
    bool PassedLambdaCuts = CheckLambdaCuts(*v0,PrimaryVertexPos,*fPIDResponse,false,RunNumber,MassInvariant,WrongMassInvariant);
    if(!PassedLambdaCuts) continue;




    float Generated_px = -999.0;
    float Generated_py = -999.0;
    float Generated_pz = -999.0;
    float Generated_px_Daughter1 = -999.0;
    float Generated_py_Daughter1 = -999.0;
    float Generated_pz_Daughter1 = -999.0;
    float Generated_px_Daughter2 = -999.0;
    float Generated_py_Daughter2 = -999.0;
    float Generated_pz_Daughter2 = -999.0;
    int PDG_Daughter1 = 0;
    int PDG_Daughter2 = 0;
    int PDG_v01 = 0;
    int PDG_v02 = 0;
    int PDG_Mother1 = 0;
    int PDG_Mother2 = 0;
    bool SameV0 = false;

    int Label_Daughter1 = TMath::Abs(AntiProtonTrack->GetLabel());
    int Label_Daughter2 = TMath::Abs(PionTrack->GetLabel());


    if(fIsMC == true){

      // get list of generated particles
      TClonesArray *MCArray = dynamic_cast <TClonesArray*> (fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!MCArray) continue;

      // Daughter1, v01, Mother1
      AliMCParticle *MCParticleDaughter1 = (AliMCParticle*) MCArray->At(Label_Daughter1);
      PDG_Daughter1 = MCParticleDaughter1->PdgCode();

      int Label_v01 = TMath::Abs(MCParticleDaughter1->GetMother());
      AliMCParticle *MCv01 = (AliMCParticle*) MCArray->At(Label_v01);
      PDG_v01 = MCv01->PdgCode();

      int Label_Mother1 = TMath::Abs(MCv01->GetMother());
      AliMCParticle *MCParticleMother1 = (AliMCParticle*) MCArray->At(Label_Mother1);
      PDG_Mother1 = MCParticleMother1->PdgCode();


      // Daughter2, v02, Mother2
      AliMCParticle *MCParticleDaughter2 = (AliMCParticle*) fMCEvent->GetTrack(Label_Daughter2);
      PDG_Daughter2 = MCParticleDaughter2->PdgCode();

      int Label_v02 = TMath::Abs(MCParticleDaughter2->GetMother());
      AliMCParticle *MCv02 = (AliMCParticle*) MCArray->At(Label_v02);
      PDG_v02 = MCv02->PdgCode();

      int Label_Mother2 = TMath::Abs(MCv02->GetMother());
      AliMCParticle *MCParticleMother2 = (AliMCParticle*) MCArray->At(Label_Mother2);
      PDG_Mother2 = MCParticleMother2->PdgCode();



      // check if daughters come from the same v0 
      SameV0 = false;
      if(Label_v01 == Label_v02) SameV0 = true;


      // where do the particles come from?
      if(MCParticleDaughter1->IsPhysicalPrimary() == true)	  PDG_v01 = 1;
      if(MCParticleDaughter1->IsSecondaryFromMaterial() == true)  PDG_v01 = 2;

      if(MCParticleDaughter2->IsPhysicalPrimary() == true)	  PDG_v02 = 1;
      if(MCParticleDaughter2->IsSecondaryFromMaterial() == true)  PDG_v02 = 2;

      if(MCv01->IsPhysicalPrimary() == true)	    PDG_Mother1 = 1;
      if(MCv01->IsSecondaryFromMaterial() == true)  PDG_Mother1 = 2;

      if(MCv02->IsPhysicalPrimary() == true)	    PDG_Mother2 = 1;
      if(MCv02->IsSecondaryFromMaterial() == true)  PDG_Mother2 = 2;

      if(SameV0 == true){

	// save generated momenta
	Generated_px = MCv01->Px();
	Generated_py = MCv01->Py();
	Generated_pz = MCv01->Pz();

      }

      Generated_px_Daughter1 = MCParticleDaughter1->Px();
      Generated_py_Daughter1 = MCParticleDaughter1->Py();
      Generated_pz_Daughter1 = MCParticleDaughter1->Pz();
      Generated_px_Daughter2 = MCParticleDaughter2->Px();
      Generated_py_Daughter2 = MCParticleDaughter2->Py();
      Generated_pz_Daughter2 = MCParticleDaughter2->Pz();


      if(DebugMC == true){

	std::cout << "\n v0: " << V0 << "/" << nV0s << std::endl;
	std::cout << "Label_Daughter1: " << Label_Daughter1 << "\t PDG_Daughter1: " << PDG_Daughter1 << "\t Label_Daughter2: " << Label_Daughter2 << "\t PDG_Daughter2: " << PDG_Daughter2 << std::endl;
	std::cout << "Label_v01: " << Label_v01 << "\t PDG_v01: " << PDG_v01 << "\t Label_v02: " << Label_v02 << "\t PDG_v02: " << PDG_v02 << "\t SameV0: " << SameV0 << std::endl;
	std::cout << "Label_Mother1: " << Label_Mother1 << "\t PDG_Mother1: " << PDG_Mother1 << "\t Label_Mother2: " << Label_Mother2 << "\t PDG_Mother2: " << PDG_Mother2 << "\n" << std::endl;

      }


    } // end of fIsMC == true







  
    TVector3 *MomentumDaughterPositive = new TVector3();
    TVector3 *MomentumDaughterNegative = new TVector3();
    TVector3 *MomentumV0 = new TVector3();
    MomentumDaughterPositive->SetXYZ(MomentumDaughterPosX,MomentumDaughterPosY,MomentumDaughterPosZ);
    MomentumDaughterNegative->SetXYZ(MomentumDaughterNegX,MomentumDaughterNegY,MomentumDaughterNegZ);
    MomentumV0->SetXYZ(MomentumV0X,MomentumV0Y,MomentumV0Z);
  
    double p_Longitudinal_Positive = MomentumDaughterPositive->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of positively charged daughter
    double p_Longitudinal_Negative = MomentumDaughterNegative->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of negatively charged daughter
  
    double Alpha = (p_Longitudinal_Positive - p_Longitudinal_Negative) / (p_Longitudinal_Positive + p_Longitudinal_Negative);
    double qT = MomentumDaughterNegative->Perp(*MomentumV0);

    

    float AntiProton_xv[2];
    float AntiProton_yv[3];
    AntiProtonTrack->GetImpactParameters(AntiProton_xv,AntiProton_yv);
    float AntiProton_DCAxy = AntiProton_xv[0];
    float AntiProton_DCAz = AntiProton_xv[1];

    double AntiProton_TOF_m2	   = -999.0;
    double AntiProton_TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus AntiProtonStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,AntiProtonTrack);
    bool AntiProtonTOFisOK = false;
    if(AntiProtonStatusTOF == AliPIDResponse::kDetPidOk) AntiProtonTOFisOK = true;
    
    if(AntiProtonTOFisOK == true){

      AntiProton_TOF_m2	    = CalculateMassSquareTOF(*AntiProtonTrack);
      AntiProton_TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(AntiProtonTrack->Pt(),AntiProton_TOF_m2,1,RunNumber);

    }

    double AntiProton_ITS_dEdx	    = -999.0;
    double AntiProton_ITS_dEdx_nSigma   = -999.0;
    int AntiProton_ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus AntiProtonStatusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,AntiProtonTrack);
    bool AntiProtonITSisOK = false;
    if(AntiProtonStatusITS == AliPIDResponse::kDetPidOk) AntiProtonITSisOK = true;
    
    if(AntiProtonITSisOK == true){

      AntiProton_ITS_dEdx	      = AntiProtonTrack->GetITSsignal();
      AntiProton_ITS_dEdx_nSigma  = CalculateSigmadEdxITS(*AntiProtonTrack,3,RunNumber);
      AntiProton_ITS_nCluster     = AntiProtonTrack->GetITSNcls();
      if(AntiProton_ITS_nCluster < 0) AntiProton_ITS_nCluster = 0;

    }



    float Pion_xv[2];
    float Pion_yv[3];
    PionTrack->GetImpactParameters(Pion_xv,Pion_yv);
    float Pion_DCAxy = Pion_xv[0];
    float Pion_DCAz = Pion_xv[1];

    double Pion_TOF_m2	   = -999.0;
    double Pion_TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus PionStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,PionTrack);
    bool PionTOFisOK = false;
    if(PionStatusTOF == AliPIDResponse::kDetPidOk) PionTOFisOK = true;
    
    if(PionTOFisOK == true){

      Pion_TOF_m2	    = CalculateMassSquareTOF(*PionTrack);
      Pion_TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(PionTrack->Pt(),Pion_TOF_m2,5,RunNumber);

    }

    double Pion_ITS_dEdx	    = -999.0;
    double Pion_ITS_dEdx_nSigma   = -999.0;
    int Pion_ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus PionStatusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,PionTrack);
    bool PionITSisOK = false;
    if(PionStatusITS == AliPIDResponse::kDetPidOk) PionITSisOK = true;
    
    if(PionITSisOK == true){

      Pion_ITS_dEdx	      = PionTrack->GetITSsignal();
      Pion_ITS_dEdx_nSigma  = CalculateSigmadEdxITS(*PionTrack,5,RunNumber);
      Pion_ITS_nCluster     = PionTrack->GetITSNcls();
      if(Pion_ITS_nCluster < 0) Pion_ITS_nCluster = 0;

    }



    AntiLambda_px		      = (float)v0->MomV0X();
    AntiLambda_py		      = (float)v0->MomV0Y();
    AntiLambda_pz		      = (float)v0->MomV0Z();
    AntiLambda_px_Generated	      = Generated_px;
    AntiLambda_py_Generated	      = Generated_py;
    AntiLambda_pz_Generated	      = Generated_pz;
    AntiLambda_Eta		      = (float)v0->PseudoRapV0();
    AntiLambda_Phi		      = (float)v0->Phi();
    AntiLambda_TransverseRadius	      = (float)v0->RadiusV0();
    AntiLambda_CosinePointingAngle    = (float)TMath::Abs(v0->CosPointingAngle(PrimaryVertexPos));
    AntiLambda_DCAv0ToPrimaryVertex   = (float)v0->DcaV0ToPrimVertex();
    AntiLambda_DCAv0Daughters	      = (float)v0->DcaV0Daughters();
    AntiLambda_Alpha		      = (float)Alpha;
    AntiLambda_qT		      = (float)qT;
    AntiLambda_DecayLength	      = (float)v0->DecayLengthV0(PrimaryVertexPos);
    AntiLambda_OpenAngle	      = (float)v0->OpenAngleV0();
    AntiLambda_PDG_Daughter1	      = PDG_Daughter1;
    AntiLambda_PDG_Daughter2	      = PDG_Daughter2;
    AntiLambda_PDG_v01		      = PDG_v01;
    AntiLambda_PDG_v02		      = PDG_v02;
    AntiLambda_PDG_Mother1	      = PDG_Mother1;
    AntiLambda_PDG_Mother2	      = PDG_Mother2;
    AntiLambda_SameV0		      = SameV0;

    AntiLambda_Daughter_AntiProton_px			    = AntiProtonTrack->Px();
    AntiLambda_Daughter_AntiProton_py			    = AntiProtonTrack->Py();
    AntiLambda_Daughter_AntiProton_pz			    = AntiProtonTrack->Pz();
    AntiLambda_Daughter_AntiProton_px_Generated		    = Generated_px_Daughter1;
    AntiLambda_Daughter_AntiProton_py_Generated		    = Generated_py_Daughter1;
    AntiLambda_Daughter_AntiProton_pz_Generated		    = Generated_pz_Daughter1;
    AntiLambda_Daughter_AntiProton_px_DecayVertex	    = MomentumDaughterNegX;
    AntiLambda_Daughter_AntiProton_py_DecayVertex	    = MomentumDaughterNegY;
    AntiLambda_Daughter_AntiProton_pz_DecayVertex	    = MomentumDaughterNegZ;
    AntiLambda_Daughter_AntiProton_pTPC			    = AntiProtonTrack->GetTPCmomentum();
    AntiLambda_Daughter_AntiProton_Eta			    = AntiProtonTrack->Eta();
    AntiLambda_Daughter_AntiProton_Phi			    = AntiProtonTrack->Phi();
    AntiLambda_Daughter_AntiProton_TPC_Chi2		    = AntiProtonTrack->GetTPCchi2();
    AntiLambda_Daughter_AntiProton_TPC_dEdx		    = AntiProtonTrack->GetTPCsignal();
    AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma	    = (float)fPIDResponse->NumberOfSigmasTPC(AntiProtonTrack,AliPID::kProton);
    AntiLambda_Daughter_AntiProton_TOF_Mass2		    = (float)AntiProton_TOF_m2;
    AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma	    = (float)AntiProton_TOF_m2_nSigma;
    AntiLambda_Daughter_AntiProton_ITS_dEdx		    = (float)AntiProton_ITS_dEdx;
    AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma	    = (float)AntiProton_ITS_dEdx_nSigma;
    AntiLambda_Daughter_AntiProton_DCAxy		    = AntiProton_DCAxy;
    AntiLambda_Daughter_AntiProton_DCAz			    = AntiProton_DCAz;
    AntiLambda_Daughter_AntiProton_TPC_nCrossedRows	    = AntiProtonTrack->GetTPCCrossedRows();
    AntiLambda_Daughter_AntiProton_TPC_nSharedCluster	    = AntiProtonTrack->GetTPCnclsS();
    AntiLambda_Daughter_AntiProton_TPC_nFindableCluster	    = AntiProtonTrack->GetTPCNclsF();
    AntiLambda_Daughter_AntiProton_TPC_nCluster		    = AntiProtonTrack->GetTPCNcls();
    AntiLambda_Daughter_AntiProton_ITS_nCluster		    = (unsigned short)AntiProton_ITS_nCluster;
    AntiLambda_Daughter_AntiProton_ID			    = Label_Daughter1;

    AntiLambda_Daughter_Pion_px			    = PionTrack->Px();
    AntiLambda_Daughter_Pion_py			    = PionTrack->Py();
    AntiLambda_Daughter_Pion_pz			    = PionTrack->Pz();
    AntiLambda_Daughter_Pion_px_Generated	    = Generated_px_Daughter2;
    AntiLambda_Daughter_Pion_py_Generated	    = Generated_py_Daughter2;
    AntiLambda_Daughter_Pion_pz_Generated	    = Generated_pz_Daughter2;
    AntiLambda_Daughter_Pion_px_DecayVertex	    = MomentumDaughterPosX;
    AntiLambda_Daughter_Pion_py_DecayVertex	    = MomentumDaughterPosY;
    AntiLambda_Daughter_Pion_pz_DecayVertex	    = MomentumDaughterPosZ;
    AntiLambda_Daughter_Pion_pTPC		    = PionTrack->GetTPCmomentum();
    AntiLambda_Daughter_Pion_Eta		    = PionTrack->Eta();
    AntiLambda_Daughter_Pion_Phi		    = PionTrack->Phi();
    AntiLambda_Daughter_Pion_TPC_Chi2		    = PionTrack->GetTPCchi2();
    AntiLambda_Daughter_Pion_TPC_dEdx		    = PionTrack->GetTPCsignal();
    AntiLambda_Daughter_Pion_TPC_dEdx_nSigma	    = (float)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
    AntiLambda_Daughter_Pion_TOF_Mass2		    = (float)Pion_TOF_m2;
    AntiLambda_Daughter_Pion_TOF_Mass2_nSigma	    = (float)Pion_TOF_m2_nSigma;
    AntiLambda_Daughter_Pion_ITS_dEdx		    = (float)Pion_ITS_dEdx;
    AntiLambda_Daughter_Pion_ITS_dEdx_nSigma	    = (float)Pion_ITS_dEdx_nSigma;
    AntiLambda_Daughter_Pion_DCAxy		    = Pion_DCAxy;
    AntiLambda_Daughter_Pion_DCAz		    = Pion_DCAz;
    AntiLambda_Daughter_Pion_TPC_nCrossedRows	    = PionTrack->GetTPCCrossedRows();
    AntiLambda_Daughter_Pion_TPC_nSharedCluster	    = PionTrack->GetTPCnclsS();
    AntiLambda_Daughter_Pion_TPC_nFindableCluster   = PionTrack->GetTPCNclsF();
    AntiLambda_Daughter_Pion_TPC_nCluster	    = PionTrack->GetTPCNcls();
    AntiLambda_Daughter_Pion_ITS_nCluster	    = (unsigned short)Pion_ITS_nCluster;
    AntiLambda_Daughter_Pion_ID			    = Label_Daughter2;


    fTempTree_AntiLambda->Fill();

    nAntiLambdasSelected++;




  } // end of AntiLambda loop




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antideuteron selection loop +++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    TTree     *fTempTree_AntiDeuteron;
    float     AntiDeuteron_px;
    float     AntiDeuteron_py;
    float     AntiDeuteron_pz;
    float     AntiDeuteron_px_Generated = 0.0;
    float     AntiDeuteron_py_Generated = 0.0;
    float     AntiDeuteron_pz_Generated = 0.0;
    float     AntiDeuteron_pTPC;
    float     AntiDeuteron_Eta;
    float     AntiDeuteron_Phi;
    float     AntiDeuteron_TPC_Chi2;
    float     AntiDeuteron_TPC_dEdx;
    float     AntiDeuteron_TPC_dEdx_nSigma;
    float     AntiDeuteron_TOF_Mass2;
    float     AntiDeuteron_TOF_Mass2_nSigma;
    float     AntiDeuteron_ITS_dEdx;
    float     AntiDeuteron_ITS_dEdx_nSigma;
    float     AntiDeuteron_DCAxy;
    float     AntiDeuteron_DCAz;
    unsigned short  AntiDeuteron_TPC_nCrossedRows;
    unsigned short  AntiDeuteron_TPC_nSharedCluster;
    unsigned short  AntiDeuteron_TPC_nFindableCluster;
    unsigned short  AntiDeuteron_TPC_nCluster;
    unsigned short  AntiDeuteron_ITS_nCluster;
    int		    AntiDeuteron_PDG;
    int		    AntiDeuteron_MotherPDG;
    unsigned int    AntiDeuteron_ID;
    unsigned long   AntiDeuteron_Event_Identifier;
    bool	    AntiDeuteron_FilterBit;

  fTempTree_AntiDeuteron = new TTree("fTempTree_AntiDeuteron","fTempTree_AntiDeuteron");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_px",&AntiDeuteron_px,"AntiDeuteron_px/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_py",&AntiDeuteron_py,"AntiDeuteron_py/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pz",&AntiDeuteron_pz,"AntiDeuteron_pz/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_px_Generated",&AntiDeuteron_px_Generated,"AntiDeuteron_px_Generated/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_py_Generated",&AntiDeuteron_py_Generated,"AntiDeuteron_py_Generated/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pz_Generated",&AntiDeuteron_pz_Generated,"AntiDeuteron_pz_Generated/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&AntiDeuteron_pTPC,"AntiDeuteron_pTPC/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&AntiDeuteron_Eta,"AntiDeuteron_Eta/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&AntiDeuteron_Phi,"AntiDeuteron_Phi/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2",&AntiDeuteron_TPC_Chi2,"AntiDeuteron_TPC_Chi2/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&AntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_nSigma",&AntiDeuteron_TPC_dEdx_nSigma,"AntiDeuteron_TPC_dEdx_nSigma/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&AntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_nSigma",&AntiDeuteron_TOF_Mass2_nSigma,"AntiDeuteron_TOF_Mass2_nSigma/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&AntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_nSigma",&AntiDeuteron_ITS_dEdx_nSigma,"AntiDeuteron_ITS_dEdx_nSigma/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&AntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&AntiDeuteron_DCAz,"AntiDeuteron_DCAz/F");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&AntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&AntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nFindableCluster",&AntiDeuteron_TPC_nFindableCluster,"AntiDeuteron_TPC_nFindableCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&AntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&AntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_PDG",&AntiDeuteron_PDG,"AntiDeuteron_PDG/I");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_MotherPDG",&AntiDeuteron_MotherPDG,"AntiDeuteron_MotherPDG/I");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_FilterBit",&AntiDeuteron_FilterBit,"AntiDeuteron_FilterBit/O");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ID",&AntiDeuteron_ID,"AntiDeuteron_ID/i");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&AntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/l");




  unsigned short nAntiDeuteronsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// antideuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antideuteron cuts
    bool PassedAntiDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,false,RunNumber);
    if(!PassedAntiDeuteronCuts) continue;




    AliMCParticle *MCParticle = 0x0;
    int Label = TMath::Abs(Track->GetLabel());

    float Generated_px = 0.0;
    float Generated_py = 0.0;
    float Generated_pz = 0.0;
    int PDG = 0;
    int MotherPDG = 0;

    if(fIsMC == true){

      MCParticle = (AliMCParticle*) fMCEvent->GetTrack(Label);
      PDG = MCParticle->PdgCode();
      if(MCParticle->IsPhysicalPrimary() == true)         MotherPDG = 1;
      if(MCParticle->IsSecondaryFromMaterial() == true)   MotherPDG = 2;
      if(MCParticle->IsSecondaryFromWeakDecay() == true){
      
        int LabelMother = TMath::Abs(MCParticle->GetMother());
        AliMCParticle *MCParticleMother = (AliMCParticle*) fMCEvent->GetTrack(LabelMother);
        MotherPDG = MCParticleMother->PdgCode();

      } 

      Generated_px = MCParticle->Px();
      Generated_py = MCParticle->Py();
      Generated_pz = MCParticle->Pz();

    } // end of fIsMC == true








  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    float DCAxy = xv[0];
    float DCAz = xv[1];

    double TOF_m2	    = -999.0;
    double TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    bool TOFisOK = false;
    if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
    if(TOFisOK){

      TOF_m2	      = CalculateMassSquareTOF(*Track);
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,4,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,4,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();
      if(ITS_nCluster < 0) ITS_nCluster = 0;

    }


    AntiDeuteron_px			= Track->Px();
    AntiDeuteron_py			= Track->Py();
    AntiDeuteron_pz			= Track->Pz();
    AntiDeuteron_px_Generated		= Generated_px;
    AntiDeuteron_py_Generated		= Generated_py;
    AntiDeuteron_pz_Generated		= Generated_pz;
    AntiDeuteron_pTPC			= Track->GetTPCmomentum();
    AntiDeuteron_Eta			= Track->Eta();
    AntiDeuteron_Phi			= Track->Phi();
    AntiDeuteron_TPC_Chi2		= Track->GetTPCchi2();
    AntiDeuteron_TPC_dEdx		= Track->GetTPCsignal();
    AntiDeuteron_TPC_dEdx_nSigma	= (float)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    AntiDeuteron_TOF_Mass2		= (float)TOF_m2;
    AntiDeuteron_TOF_Mass2_nSigma	= (float)TOF_m2_nSigma;
    AntiDeuteron_ITS_dEdx		= (float)ITS_dEdx;
    AntiDeuteron_ITS_dEdx_nSigma	= (float)ITS_dEdx_nSigma;
    AntiDeuteron_DCAxy			= DCAxy;
    AntiDeuteron_DCAz			= DCAz;
    AntiDeuteron_TPC_nCrossedRows	= Track->GetTPCCrossedRows();
    AntiDeuteron_TPC_nSharedCluster	= Track->GetTPCnclsS();
    AntiDeuteron_TPC_nFindableCluster	= Track->GetTPCNclsF();
    AntiDeuteron_TPC_nCluster		= Track->GetTPCNcls();
    AntiDeuteron_ITS_nCluster		= (unsigned short)ITS_nCluster;
    AntiDeuteron_ID			= track;
    AntiDeuteron_Event_Identifier	= EventID;
    AntiDeuteron_PDG			= PDG;
    AntiDeuteron_MotherPDG		= MotherPDG;
    AntiDeuteron_FilterBit		= Track->TestFilterBit(BIT(0));
 
    fTempTree_AntiDeuteron->Fill();
    nAntiDeuteronsSelected++;

  } // end of antideuteron loop




  if((fSavePairsOnly == false) || ((fSavePairsOnly == true) && (nAntiLambdasSelected > 0) && (nAntiDeuteronsSelected > 0))){

    for(int AntiLambda = 0; AntiLambda < nAntiLambdasSelected; AntiLambda++){


      TBranch *Branch_AntiLambda_px			  = fTempTree_AntiLambda->GetBranch("AntiLambda_px");
      TBranch *Branch_AntiLambda_py			  = fTempTree_AntiLambda->GetBranch("AntiLambda_py");
      TBranch *Branch_AntiLambda_pz			  = fTempTree_AntiLambda->GetBranch("AntiLambda_pz");
      TBranch *Branch_AntiLambda_px_Generated		  = fTempTree_AntiLambda->GetBranch("AntiLambda_px_Generated");
      TBranch *Branch_AntiLambda_py_Generated		  = fTempTree_AntiLambda->GetBranch("AntiLambda_py_Generated");
      TBranch *Branch_AntiLambda_pz_Generated		  = fTempTree_AntiLambda->GetBranch("AntiLambda_pz_Generated");
      TBranch *Branch_AntiLambda_Eta			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Eta");
      TBranch *Branch_AntiLambda_Phi			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Phi");
      TBranch *Branch_AntiLambda_TransverseRadius	  = fTempTree_AntiLambda->GetBranch("AntiLambda_TransverseRadius");
      TBranch *Branch_AntiLambda_CosinePointingAngle	  = fTempTree_AntiLambda->GetBranch("AntiLambda_CosinePointingAngle");
      TBranch *Branch_AntiLambda_DCAv0ToPrimaryVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_DCAv0ToPrimaryVertex");
      TBranch *Branch_AntiLambda_DCAv0Daughters		  = fTempTree_AntiLambda->GetBranch("AntiLambda_DCAv0Daughters");
      TBranch *Branch_AntiLambda_Alpha			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Alpha");
      TBranch *Branch_AntiLambda_qT			  = fTempTree_AntiLambda->GetBranch("AntiLambda_qT");
      TBranch *Branch_AntiLambda_DecayLength		  = fTempTree_AntiLambda->GetBranch("AntiLambda_DecayLength");
      TBranch *Branch_AntiLambda_OpenAngle		  = fTempTree_AntiLambda->GetBranch("AntiLambda_OpenAngle");
      TBranch *Branch_AntiLambda_PDG_Daughter1		  = fTempTree_AntiLambda->GetBranch("AntiLambda_PDG_Daughter1");
      TBranch *Branch_AntiLambda_PDG_Daughter2		  = fTempTree_AntiLambda->GetBranch("AntiLambda_PDG_Daughter2");
      TBranch *Branch_AntiLambda_PDG_v01		  = fTempTree_AntiLambda->GetBranch("AntiLambda_PDG_v01");
      TBranch *Branch_AntiLambda_PDG_v02		  = fTempTree_AntiLambda->GetBranch("AntiLambda_PDG_v02");
      TBranch *Branch_AntiLambda_PDG_Mother1		  = fTempTree_AntiLambda->GetBranch("AntiLambda_PDG_Mother1");
      TBranch *Branch_AntiLambda_PDG_Mother2		  = fTempTree_AntiLambda->GetBranch("AntiLambda_PDG_Mother2");
      TBranch *Branch_AntiLambda_SameV0			  = fTempTree_AntiLambda->GetBranch("AntiLambda_SameV0");

      TBranch *Branch_AntiLambda_Daughter_AntiProton_px			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_px");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_py			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_py");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_pz			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_pz");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_px_Generated	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_px_Generated");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_py_Generated	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_py_Generated");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_pz_Generated	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_pz_Generated");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_px_DecayVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_px_DecayVertex");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_py_DecayVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_py_DecayVertex");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_pz_DecayVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_pz_DecayVertex");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_pTPC		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_pTPC");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_Eta		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_Eta");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_Phi		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_Phi");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_Chi2		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_Chi2");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_dEdx		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_dEdx");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TOF_Mass2		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TOF_Mass2");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_ITS_dEdx		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_ITS_dEdx");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_DCAxy		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_DCAxy");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_DCAz		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_DCAz");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_nCrossedRows	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_nCrossedRows");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_nSharedCluster	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_nSharedCluster");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_nFindableCluster = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_nFindableCluster");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_TPC_nCluster	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_TPC_nCluster");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_ITS_nCluster	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_ITS_nCluster");
      TBranch *Branch_AntiLambda_Daughter_AntiProton_ID			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_AntiProton_ID");

      TBranch *Branch_AntiLambda_Daughter_Pion_px		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_px");
      TBranch *Branch_AntiLambda_Daughter_Pion_py		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_py");
      TBranch *Branch_AntiLambda_Daughter_Pion_pz		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_pz");
      TBranch *Branch_AntiLambda_Daughter_Pion_px_Generated	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_px_Generated");
      TBranch *Branch_AntiLambda_Daughter_Pion_py_Generated	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_py_Generated");
      TBranch *Branch_AntiLambda_Daughter_Pion_pz_Generated	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_pz_Generated");
      TBranch *Branch_AntiLambda_Daughter_Pion_px_DecayVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_px_DecayVertex");
      TBranch *Branch_AntiLambda_Daughter_Pion_py_DecayVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_py_DecayVertex");
      TBranch *Branch_AntiLambda_Daughter_Pion_pz_DecayVertex	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_pz_DecayVertex");
      TBranch *Branch_AntiLambda_Daughter_Pion_pTPC		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_pTPC");
      TBranch *Branch_AntiLambda_Daughter_Pion_Eta		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_Eta");
      TBranch *Branch_AntiLambda_Daughter_Pion_Phi		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_Phi");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_Chi2		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_Chi2");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_dEdx		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_dEdx");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_dEdx_nSigma	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_dEdx_nSigma");
      TBranch *Branch_AntiLambda_Daughter_Pion_TOF_Mass2		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TOF_Mass2");
      TBranch *Branch_AntiLambda_Daughter_Pion_TOF_Mass2_nSigma	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TOF_Mass2_nSigma");
      TBranch *Branch_AntiLambda_Daughter_Pion_ITS_dEdx		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_ITS_dEdx");
      TBranch *Branch_AntiLambda_Daughter_Pion_ITS_dEdx_nSigma	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_ITS_dEdx_nSigma");
      TBranch *Branch_AntiLambda_Daughter_Pion_DCAxy		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_DCAxy");
      TBranch *Branch_AntiLambda_Daughter_Pion_DCAz		  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_DCAz");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_nCrossedRows	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_nCrossedRows");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_nSharedCluster	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_nSharedCluster");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_nFindableCluster = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_nFindableCluster");
      TBranch *Branch_AntiLambda_Daughter_Pion_TPC_nCluster	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_TPC_nCluster");
      TBranch *Branch_AntiLambda_Daughter_Pion_ITS_nCluster	  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_ITS_nCluster");
      TBranch *Branch_AntiLambda_Daughter_Pion_ID			  = fTempTree_AntiLambda->GetBranch("AntiLambda_Daughter_Pion_ID");


      Branch_AntiLambda_px->SetAddress(&fAntiLambda_px);
      Branch_AntiLambda_py->SetAddress(&fAntiLambda_py);
      Branch_AntiLambda_pz->SetAddress(&fAntiLambda_pz);
      Branch_AntiLambda_px_Generated->SetAddress(&fAntiLambda_px_Generated);
      Branch_AntiLambda_py_Generated->SetAddress(&fAntiLambda_py_Generated);
      Branch_AntiLambda_pz_Generated->SetAddress(&fAntiLambda_pz_Generated);
      Branch_AntiLambda_Eta->SetAddress(&fAntiLambda_Eta);
      Branch_AntiLambda_Phi->SetAddress(&fAntiLambda_Phi);
      Branch_AntiLambda_TransverseRadius->SetAddress(&fAntiLambda_TransverseRadius);
      Branch_AntiLambda_CosinePointingAngle->SetAddress(&fAntiLambda_CosinePointingAngle);
      Branch_AntiLambda_DCAv0ToPrimaryVertex->SetAddress(&fAntiLambda_DCAv0ToPrimaryVertex);
      Branch_AntiLambda_DCAv0Daughters->SetAddress(&fAntiLambda_DCAv0Daughters);
      Branch_AntiLambda_Alpha->SetAddress(&fAntiLambda_Alpha);
      Branch_AntiLambda_qT->SetAddress(&fAntiLambda_qT);
      Branch_AntiLambda_DecayLength->SetAddress(&fAntiLambda_DecayLength);
      Branch_AntiLambda_OpenAngle->SetAddress(&fAntiLambda_OpenAngle);
      Branch_AntiLambda_PDG_Daughter1->SetAddress(&fAntiLambda_PDG_Daughter1);
      Branch_AntiLambda_PDG_Daughter2->SetAddress(&fAntiLambda_PDG_Daughter2);
      Branch_AntiLambda_PDG_v01->SetAddress(&fAntiLambda_PDG_v01);
      Branch_AntiLambda_PDG_v02->SetAddress(&fAntiLambda_PDG_v02);
      Branch_AntiLambda_PDG_Mother1->SetAddress(&fAntiLambda_PDG_Mother1);
      Branch_AntiLambda_PDG_Mother2->SetAddress(&fAntiLambda_PDG_Mother2);
      Branch_AntiLambda_SameV0->SetAddress(&fAntiLambda_SameV0);

      Branch_AntiLambda_Daughter_AntiProton_px->SetAddress(&fAntiLambda_Daughter_AntiProton_px);
      Branch_AntiLambda_Daughter_AntiProton_py->SetAddress(&fAntiLambda_Daughter_AntiProton_py);
      Branch_AntiLambda_Daughter_AntiProton_pz->SetAddress(&fAntiLambda_Daughter_AntiProton_pz);
      Branch_AntiLambda_Daughter_AntiProton_px_Generated->SetAddress(&fAntiLambda_Daughter_AntiProton_px_Generated);
      Branch_AntiLambda_Daughter_AntiProton_py_Generated->SetAddress(&fAntiLambda_Daughter_AntiProton_py_Generated);
      Branch_AntiLambda_Daughter_AntiProton_pz_Generated->SetAddress(&fAntiLambda_Daughter_AntiProton_pz_Generated);
      Branch_AntiLambda_Daughter_AntiProton_px_DecayVertex->SetAddress(&fAntiLambda_Daughter_AntiProton_px_DecayVertex);
      Branch_AntiLambda_Daughter_AntiProton_py_DecayVertex->SetAddress(&fAntiLambda_Daughter_AntiProton_py_DecayVertex);
      Branch_AntiLambda_Daughter_AntiProton_pz_DecayVertex->SetAddress(&fAntiLambda_Daughter_AntiProton_pz_DecayVertex);
      Branch_AntiLambda_Daughter_AntiProton_pTPC->SetAddress(&fAntiLambda_Daughter_AntiProton_pTPC);
      Branch_AntiLambda_Daughter_AntiProton_Eta->SetAddress(&fAntiLambda_Daughter_AntiProton_Eta);
      Branch_AntiLambda_Daughter_AntiProton_Phi->SetAddress(&fAntiLambda_Daughter_AntiProton_Phi);
      Branch_AntiLambda_Daughter_AntiProton_TPC_Chi2->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_Chi2);
      Branch_AntiLambda_Daughter_AntiProton_TPC_dEdx->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_dEdx);
      Branch_AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma);
      Branch_AntiLambda_Daughter_AntiProton_TOF_Mass2->SetAddress(&fAntiLambda_Daughter_AntiProton_TOF_Mass2);
      Branch_AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma->SetAddress(&fAntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma);
      Branch_AntiLambda_Daughter_AntiProton_ITS_dEdx->SetAddress(&fAntiLambda_Daughter_AntiProton_ITS_dEdx);
      Branch_AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma->SetAddress(&fAntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma);
      Branch_AntiLambda_Daughter_AntiProton_DCAxy->SetAddress(&fAntiLambda_Daughter_AntiProton_DCAxy);
      Branch_AntiLambda_Daughter_AntiProton_DCAz->SetAddress(&fAntiLambda_Daughter_AntiProton_DCAz);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nCrossedRows->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nSharedCluster->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_nSharedCluster);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nFindableCluster->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nCluster->SetAddress(&fAntiLambda_Daughter_AntiProton_TPC_nCluster);
      Branch_AntiLambda_Daughter_AntiProton_ITS_nCluster->SetAddress(&fAntiLambda_Daughter_AntiProton_ITS_nCluster);
      Branch_AntiLambda_Daughter_AntiProton_ID->SetAddress(&fAntiLambda_Daughter_AntiProton_ID);

      Branch_AntiLambda_Daughter_Pion_px->SetAddress(&fAntiLambda_Daughter_Pion_px);
      Branch_AntiLambda_Daughter_Pion_py->SetAddress(&fAntiLambda_Daughter_Pion_py);
      Branch_AntiLambda_Daughter_Pion_pz->SetAddress(&fAntiLambda_Daughter_Pion_pz);
      Branch_AntiLambda_Daughter_Pion_px_Generated->SetAddress(&fAntiLambda_Daughter_Pion_px_Generated);
      Branch_AntiLambda_Daughter_Pion_py_Generated->SetAddress(&fAntiLambda_Daughter_Pion_py_Generated);
      Branch_AntiLambda_Daughter_Pion_pz_Generated->SetAddress(&fAntiLambda_Daughter_Pion_pz_Generated);
      Branch_AntiLambda_Daughter_Pion_px_DecayVertex->SetAddress(&fAntiLambda_Daughter_Pion_px_DecayVertex);
      Branch_AntiLambda_Daughter_Pion_py_DecayVertex->SetAddress(&fAntiLambda_Daughter_Pion_py_DecayVertex);
      Branch_AntiLambda_Daughter_Pion_pz_DecayVertex->SetAddress(&fAntiLambda_Daughter_Pion_pz_DecayVertex);
      Branch_AntiLambda_Daughter_Pion_pTPC->SetAddress(&fAntiLambda_Daughter_Pion_pTPC);
      Branch_AntiLambda_Daughter_Pion_Eta->SetAddress(&fAntiLambda_Daughter_Pion_Eta);
      Branch_AntiLambda_Daughter_Pion_Phi->SetAddress(&fAntiLambda_Daughter_Pion_Phi);
      Branch_AntiLambda_Daughter_Pion_TPC_Chi2->SetAddress(&fAntiLambda_Daughter_Pion_TPC_Chi2);
      Branch_AntiLambda_Daughter_Pion_TPC_dEdx->SetAddress(&fAntiLambda_Daughter_Pion_TPC_dEdx);
      Branch_AntiLambda_Daughter_Pion_TPC_dEdx_nSigma->SetAddress(&fAntiLambda_Daughter_Pion_TPC_dEdx_nSigma);
      Branch_AntiLambda_Daughter_Pion_TOF_Mass2->SetAddress(&fAntiLambda_Daughter_Pion_TOF_Mass2);
      Branch_AntiLambda_Daughter_Pion_TOF_Mass2_nSigma->SetAddress(&fAntiLambda_Daughter_Pion_TOF_Mass2_nSigma);
      Branch_AntiLambda_Daughter_Pion_ITS_dEdx->SetAddress(&fAntiLambda_Daughter_Pion_ITS_dEdx);
      Branch_AntiLambda_Daughter_Pion_ITS_dEdx_nSigma->SetAddress(&fAntiLambda_Daughter_Pion_ITS_dEdx_nSigma);
      Branch_AntiLambda_Daughter_Pion_DCAxy->SetAddress(&fAntiLambda_Daughter_Pion_DCAxy);
      Branch_AntiLambda_Daughter_Pion_DCAz->SetAddress(&fAntiLambda_Daughter_Pion_DCAz);
      Branch_AntiLambda_Daughter_Pion_TPC_nCrossedRows->SetAddress(&fAntiLambda_Daughter_Pion_TPC_nCrossedRows);
      Branch_AntiLambda_Daughter_Pion_TPC_nSharedCluster->SetAddress(&fAntiLambda_Daughter_Pion_TPC_nSharedCluster);
      Branch_AntiLambda_Daughter_Pion_TPC_nFindableCluster->SetAddress(&fAntiLambda_Daughter_Pion_TPC_nFindableCluster);
      Branch_AntiLambda_Daughter_Pion_TPC_nCluster->SetAddress(&fAntiLambda_Daughter_Pion_TPC_nCluster);
      Branch_AntiLambda_Daughter_Pion_ITS_nCluster->SetAddress(&fAntiLambda_Daughter_Pion_ITS_nCluster);
      Branch_AntiLambda_Daughter_Pion_ID->SetAddress(&fAntiLambda_Daughter_Pion_ID);



      Branch_AntiLambda_px->SetAutoDelete(true);
      Branch_AntiLambda_py->SetAutoDelete(true);
      Branch_AntiLambda_pz->SetAutoDelete(true);
      Branch_AntiLambda_px_Generated->SetAutoDelete(true);
      Branch_AntiLambda_py_Generated->SetAutoDelete(true);
      Branch_AntiLambda_pz_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Eta->SetAutoDelete(true);
      Branch_AntiLambda_Phi->SetAutoDelete(true);
      Branch_AntiLambda_TransverseRadius->SetAutoDelete(true);
      Branch_AntiLambda_CosinePointingAngle->SetAutoDelete(true);
      Branch_AntiLambda_DCAv0ToPrimaryVertex->SetAutoDelete(true);
      Branch_AntiLambda_DCAv0Daughters->SetAutoDelete(true);
      Branch_AntiLambda_Alpha->SetAutoDelete(true);
      Branch_AntiLambda_qT->SetAutoDelete(true);
      Branch_AntiLambda_DecayLength->SetAutoDelete(true);
      Branch_AntiLambda_OpenAngle->SetAutoDelete(true);
      Branch_AntiLambda_PDG_Daughter1->SetAutoDelete(true);
      Branch_AntiLambda_PDG_Daughter2->SetAutoDelete(true);
      Branch_AntiLambda_PDG_v01->SetAutoDelete(true);
      Branch_AntiLambda_PDG_v02->SetAutoDelete(true);
      Branch_AntiLambda_PDG_Mother1->SetAutoDelete(true);
      Branch_AntiLambda_PDG_Mother2->SetAutoDelete(true);
      Branch_AntiLambda_SameV0->SetAutoDelete(true);

      Branch_AntiLambda_Daughter_AntiProton_px->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_py->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_pz->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_px_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_py_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_pz_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_px_DecayVertex->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_py_DecayVertex->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_pz_DecayVertex->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_pTPC->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_Eta->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_Phi->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_Chi2->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TOF_Mass2->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_ITS_dEdx->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_DCAxy->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_DCAz->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nFindableCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_TPC_nCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_ITS_nCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_AntiProton_ID->SetAutoDelete(true);

      Branch_AntiLambda_Daughter_Pion_px->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_py->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_pz->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_px_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_py_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_pz_Generated->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_px_DecayVertex->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_py_DecayVertex->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_pz_DecayVertex->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_pTPC->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_Eta->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_Phi->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_Chi2->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TOF_Mass2->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_ITS_dEdx->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_DCAxy->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_DCAz->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_nFindableCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_TPC_nCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_ITS_nCluster->SetAutoDelete(true);
      Branch_AntiLambda_Daughter_Pion_ID->SetAutoDelete(true);

      fTempTree_AntiLambda->GetEntry(AntiLambda);


      fAntiLambda_Event_Multiplicity    = Multiplicity;
      fAntiLambda_Event_Centrality	= Centrality;
      fAntiLambda_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fAntiLambda_Event_BField		= BField;
      fAntiLambda_Event_Identifier	= EventID;
      fAntiLambda_Event_ParticleNumber	= AntiLambda;
      fAntiLambda_Event_nParticles	= nAntiLambdasSelected;
      

      fSaveTree_AntiLambda->Fill();


    } // end of loop (copy AntiLambda)


    for(int AntiDeuteron = 0; AntiDeuteron < nAntiDeuteronsSelected; AntiDeuteron++){

      TBranch *Branch_AntiDeuteron_px			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_px");
      TBranch *Branch_AntiDeuteron_py			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_py");
      TBranch *Branch_AntiDeuteron_pz			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pz");
      TBranch *Branch_AntiDeuteron_px_Generated		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_px_Generated");
      TBranch *Branch_AntiDeuteron_py_Generated		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_py_Generated");
      TBranch *Branch_AntiDeuteron_pz_Generated		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pz_Generated");
      TBranch *Branch_AntiDeuteron_pTPC			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pTPC");
      TBranch *Branch_AntiDeuteron_Eta			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Eta");
      TBranch *Branch_AntiDeuteron_Phi			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Phi");
      TBranch *Branch_AntiDeuteron_TPC_Chi2		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_Chi2");
      TBranch *Branch_AntiDeuteron_TPC_dEdx		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_dEdx");
      TBranch *Branch_AntiDeuteron_TPC_dEdx_nSigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_dEdx_nSigma");
      TBranch *Branch_AntiDeuteron_TOF_Mass2		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Mass2");
      TBranch *Branch_AntiDeuteron_TOF_Mass2_nSigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Mass2_nSigma");
      TBranch *Branch_AntiDeuteron_ITS_dEdx		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ITS_dEdx");
      TBranch *Branch_AntiDeuteron_ITS_dEdx_nSigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ITS_dEdx_nSigma");
      TBranch *Branch_AntiDeuteron_DCAxy		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_DCAxy");
      TBranch *Branch_AntiDeuteron_DCAz			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_DCAz");
      TBranch *Branch_AntiDeuteron_TPC_nCrossedRows	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nCrossedRows");
      TBranch *Branch_AntiDeuteron_TPC_nSharedCluster   = fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nSharedCluster");
      TBranch *Branch_AntiDeuteron_TPC_nFindableCluster = fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nFindableCluster");
      TBranch *Branch_AntiDeuteron_TPC_nCluster		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nCluster");
      TBranch *Branch_AntiDeuteron_ITS_nCluster		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ITS_nCluster");
      TBranch *Branch_AntiDeuteron_PDG			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_PDG");
      TBranch *Branch_AntiDeuteron_MotherPDG		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_MotherPDG");
      TBranch *Branch_AntiDeuteron_FilterBit		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_FilterBit");
      TBranch *Branch_AntiDeuteron_ID			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ID");
      TBranch *Branch_AntiDeuteron_Event_Identifier	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Event_Identifier");


      Branch_AntiDeuteron_px->SetAddress(&fAntiDeuteron_px);
      Branch_AntiDeuteron_py->SetAddress(&fAntiDeuteron_py);
      Branch_AntiDeuteron_pz->SetAddress(&fAntiDeuteron_pz);
      Branch_AntiDeuteron_px_Generated->SetAddress(&fAntiDeuteron_px_Generated);
      Branch_AntiDeuteron_py_Generated->SetAddress(&fAntiDeuteron_py_Generated);
      Branch_AntiDeuteron_pz_Generated->SetAddress(&fAntiDeuteron_pz_Generated);
      Branch_AntiDeuteron_pTPC->SetAddress(&fAntiDeuteron_pTPC);
      Branch_AntiDeuteron_Eta->SetAddress(&fAntiDeuteron_Eta);
      Branch_AntiDeuteron_Phi->SetAddress(&fAntiDeuteron_Phi);
      Branch_AntiDeuteron_TPC_Chi2->SetAddress(&fAntiDeuteron_TPC_Chi2);
      Branch_AntiDeuteron_TPC_dEdx->SetAddress(&fAntiDeuteron_TPC_dEdx);
      Branch_AntiDeuteron_TPC_dEdx_nSigma->SetAddress(&fAntiDeuteron_TPC_dEdx_nSigma);
      Branch_AntiDeuteron_TOF_Mass2->SetAddress(&fAntiDeuteron_TOF_Mass2);
      Branch_AntiDeuteron_TOF_Mass2_nSigma->SetAddress(&fAntiDeuteron_TOF_Mass2_nSigma);
      Branch_AntiDeuteron_ITS_dEdx->SetAddress(&fAntiDeuteron_ITS_dEdx);
      Branch_AntiDeuteron_ITS_dEdx_nSigma->SetAddress(&fAntiDeuteron_ITS_dEdx_nSigma);
      Branch_AntiDeuteron_DCAxy->SetAddress(&fAntiDeuteron_DCAxy);
      Branch_AntiDeuteron_DCAz->SetAddress(&fAntiDeuteron_DCAz);
      Branch_AntiDeuteron_TPC_nCrossedRows->SetAddress(&fAntiDeuteron_TPC_nCrossedRows);
      Branch_AntiDeuteron_TPC_nSharedCluster->SetAddress(&fAntiDeuteron_TPC_nSharedCluster);
      Branch_AntiDeuteron_TPC_nFindableCluster->SetAddress(&fAntiDeuteron_TPC_nFindableCluster);
      Branch_AntiDeuteron_TPC_nCluster->SetAddress(&fAntiDeuteron_TPC_nCluster);
      Branch_AntiDeuteron_ITS_nCluster->SetAddress(&fAntiDeuteron_ITS_nCluster);
      Branch_AntiDeuteron_PDG->SetAddress(&fAntiDeuteron_PDG);
      Branch_AntiDeuteron_MotherPDG->SetAddress(&fAntiDeuteron_MotherPDG);
      Branch_AntiDeuteron_FilterBit->SetAddress(&fAntiDeuteron_FilterBit);
      Branch_AntiDeuteron_ID->SetAddress(&fAntiDeuteron_ID);
      Branch_AntiDeuteron_Event_Identifier->SetAddress(&fAntiDeuteron_Event_Identifier);


      Branch_AntiDeuteron_px->SetAutoDelete(true);
      Branch_AntiDeuteron_py->SetAutoDelete(true);
      Branch_AntiDeuteron_pz->SetAutoDelete(true);
      Branch_AntiDeuteron_px_Generated->SetAutoDelete(true);
      Branch_AntiDeuteron_py_Generated->SetAutoDelete(true);
      Branch_AntiDeuteron_pz_Generated->SetAutoDelete(true);
      Branch_AntiDeuteron_pTPC->SetAutoDelete(true);
      Branch_AntiDeuteron_Eta->SetAutoDelete(true);
      Branch_AntiDeuteron_Phi->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_Chi2->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Mass2->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_AntiDeuteron_ITS_dEdx->SetAutoDelete(true);
      Branch_AntiDeuteron_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiDeuteron_DCAxy->SetAutoDelete(true);
      Branch_AntiDeuteron_DCAz->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nFindableCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_ITS_nCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_PDG->SetAutoDelete(true);
      Branch_AntiDeuteron_MotherPDG->SetAutoDelete(true);
      Branch_AntiDeuteron_FilterBit->SetAutoDelete(true);
      Branch_AntiDeuteron_ID->SetAutoDelete(true);
      Branch_AntiDeuteron_Event_Identifier->SetAutoDelete(true);

      fTempTree_AntiDeuteron->GetEntry(AntiDeuteron);


      fAntiDeuteron_Event_Multiplicity	  = Multiplicity;
      fAntiDeuteron_Event_Centrality	  = Centrality;
      fAntiDeuteron_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fAntiDeuteron_Event_BField	  = BField;
      fAntiDeuteron_Event_ParticleNumber  = AntiDeuteron;
      fAntiDeuteron_Event_nParticles	  = nAntiDeuteronsSelected;


      fSaveTree_AntiDeuteron->Fill();


    } // end of loop (copy antideuterons)


  } // end of same-event

  fTempTree_AntiLambda->Delete();
  fTempTree_AntiDeuteron->Delete();












  PostData(1,fSaveTree_Lambda);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiLambda);
  PostData(4,fSaveTree_AntiDeuteron);
  PostData(5,fHistoList);

} // end of UserExec















void AliAnalysisTask_Ld_CreateTrees_PairsOnly::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
double AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateBetaTOF(AliAODTrack &track)
{

  double length = track.GetIntegratedLength(); // cm
  if(TMath::IsNaN(length)) return -999.0;
  if(length <= 350.0) return -999.0;

  double c = TMath::C(); // m/s
  double end_time = track.GetTOFsignal(); // ps
  double start_time = fPIDResponse->GetTOFResponse().GetStartTime(track.GetTPCmomentum()); // ps

  if(TMath::IsNaN(end_time)) return -999.0;
  if(TMath::IsNaN(start_time)) return -999.0;

  double time = (end_time - start_time) * 1e-12; // ps -> s
  double velocity = (length*0.01) / time; // m/s
  double beta = velocity / c;

  return beta;

} // end of CalculateBetaTOF






// calculate the mass² TOF
double AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateMassSquareTOF(AliAODTrack &track)
{

  double mass2 = -999.0;

  double p = track.P();
  double beta = CalculateBetaTOF(track);

  if(TMath::IsNaN(p)) return mass2;
  if(TMath::IsNaN(beta)) return mass2;

  if(beta > 0.0){

    mass2 = (1/(beta*beta)-1) * (p*p);

  }

  return mass2;

} // end of CalculateMassSquareTOF








double AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateSigmaMassSquareTOF(double pT, double massSq, int ParticleSpecies, int RunNumber)
{

  double SigmaParticle = -999.0;
  if(massSq < -990.0) return SigmaParticle;


  bool MetaLHC16 = false;
  bool MetaLHC17 = false;
  bool MetaLHC18 = false;
  bool LHC18q = false;
  bool LHC18r = false;

  if((RunNumber >= 252235) && (RunNumber <= 264347)) MetaLHC16 = true;
  if((RunNumber >= 270581) && (RunNumber <= 282704)) MetaLHC17 = true;
  if((RunNumber >= 285009) && (RunNumber <= 294925)) MetaLHC18 = true;
  if((RunNumber >= 295585) && (RunNumber <= 296623)) LHC18q = true;
  if((RunNumber >= 296690) && (RunNumber <= 297585)) LHC18r = true;


  bool isProton	      = false;
  bool isDeuteron     = false;
  bool isAntiProton   = false;
  bool isAntiDeuteron = false;
  bool isPion	      = false;
  bool isAntiPion     = false;

  if(ParticleSpecies == 1) isProton	  = true;
  if(ParticleSpecies == 2) isDeuteron	  = true;
  if(ParticleSpecies == 3) isAntiProton	  = true;
  if(ParticleSpecies == 4) isAntiDeuteron = true;
  if(ParticleSpecies == 5) isPion	  = true;
  if(ParticleSpecies == 6) isAntiPion	  = true;
  

  TF1 *Mean = new TF1("Mean","[0] + ([1] * (x)) + [2] * pow((1 -([3] / (x))),[4])",0.0,6.0);
  TF1 *Sigma = new TF1("Sigma","[0] + ([1] * (x)) + [2] * pow((1 -([3] / (x))),[4])",0.0,6.0);

  // MetaLHC16 pp data 
  if((MetaLHC16 == true) && (isProton == true)){

    Mean->FixParameter(0,0.886844);
    Mean->FixParameter(1,0.00635951);
    Mean->FixParameter(2,-1.70888e-06);
    Mean->FixParameter(3,22.3206);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,383.247);
    Sigma->FixParameter(1,0.0635129);
    Sigma->FixParameter(2,-383.326);
    Sigma->FixParameter(3,0.00118195);
    Sigma->FixParameter(4,0.108687);

  }


  if((MetaLHC16 == true) && (isDeuteron == true)){

    Mean->FixParameter(0,3.53557);
    Mean->FixParameter(1,0.00503793);
    Mean->FixParameter(2,-0.000335782);
    Mean->FixParameter(3,10.3859);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,1.02977);
    Sigma->FixParameter(1,0.153259);
    Sigma->FixParameter(2,-1.41417);
    Sigma->FixParameter(3,0.0996178);
    Sigma->FixParameter(4,3);

  }

  if((MetaLHC16 == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.885019);
    Mean->FixParameter(1,0.00671772);
    Mean->FixParameter(2,-1.45866e-06);
    Mean->FixParameter(3,23.4559);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,498.345);
    Sigma->FixParameter(1,0.0771085);
    Sigma->FixParameter(2,-498.456);
    Sigma->FixParameter(3,0.00100001);
    Sigma->FixParameter(4,0.13076);

  }

  if((MetaLHC16 == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,3.53154);
    Mean->FixParameter(1,0.0050005);
    Mean->FixParameter(2,-0.000985306);
    Mean->FixParameter(3,7.57384);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.870801);
    Sigma->FixParameter(1,0.088904);
    Sigma->FixParameter(2,-1.00007);
    Sigma->FixParameter(3,0.0530252);
    Sigma->FixParameter(4,3);

  }

  if((MetaLHC16 == true) && (isPion == true)){

    Mean->FixParameter(0,0.0155002);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.00550021);
    Mean->FixParameter(3,-3.22578e+08);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,-0.00402705);
    Sigma->FixParameter(1,0.0177918);
    Sigma->FixParameter(2,-0.00101698);
    Sigma->FixParameter(3,-1.89407);
    Sigma->FixParameter(4,-27.8531);

  }

  if((MetaLHC16 == true) && (isAntiPion == true)){

    Mean->FixParameter(0,0.0182403);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.0032403);
    Mean->FixParameter(3,-1.59699e+07);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,-0.0684207);
    Sigma->FixParameter(1,0.0168145);
    Sigma->FixParameter(2,0.0621792);
    Sigma->FixParameter(3,0.000471401);
    Sigma->FixParameter(4,-40.3507);

  }


  // MetaLHC17 pp data 
  if((MetaLHC17 == true) && (isProton == true)){

    Mean->FixParameter(0,0.889429);
    Mean->FixParameter(1,0.00500005);
    Mean->FixParameter(2,-5.08905e-07);
    Mean->FixParameter(3,33.165);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.0171245);
    Sigma->FixParameter(1,0.0697961);
    Sigma->FixParameter(2,-0.168793);
    Sigma->FixParameter(3,0.0270737);
    Sigma->FixParameter(4,39.7396);

  }

  if((MetaLHC17 == true) && (isDeuteron == true)){

    Mean->FixParameter(0,3.55905);
    Mean->FixParameter(1,-0.000832746);
    Mean->FixParameter(2,-1.64011e-05);
    Mean->FixParameter(3,25.9856);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,1.92524);
    Sigma->FixParameter(1,0.0938953);
    Sigma->FixParameter(2,-2.12169);
    Sigma->FixParameter(3,0.0355014);
    Sigma->FixParameter(4,3);

  }

  if((MetaLHC17 == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.885178);
    Mean->FixParameter(1,0.00500141);
    Mean->FixParameter(2,-7.18764e-06);
    Mean->FixParameter(3,14.0775);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.00380517);
    Sigma->FixParameter(1,0.0767935);
    Sigma->FixParameter(2,-0.19083);
    Sigma->FixParameter(3,0.0193695);
    Sigma->FixParameter(4,65.9204);

  }


  if((MetaLHC17 == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,3.61532);
    Mean->FixParameter(1,-0.0254348);
    Mean->FixParameter(2,-0.0191476);
    Mean->FixParameter(3,3.31381);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.0113851);
    Sigma->FixParameter(1,0.0442788);
    Sigma->FixParameter(2,-0.00524713);
    Sigma->FixParameter(3,3.4115);
    Sigma->FixParameter(4,3);

  }

  if((MetaLHC17 == true) && (isPion == true)){

    Mean->FixParameter(0,0.0155);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.0055);
    Mean->FixParameter(3,-3e+09);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }

  if((MetaLHC17 == true) && (isAntiPion == true)){

    Mean->FixParameter(0,0.0172665);
    Mean->FixParameter(1,0.00177955);
    Mean->FixParameter(2,0.00226655);
    Mean->FixParameter(3,-3.03171e+07);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }

  // MetaLHC18 pp data 
  if((MetaLHC18 == true) && (isProton == true)){

    Mean->FixParameter(0,0.887931);
    Mean->FixParameter(1,0.00574878);
    Mean->FixParameter(2,-8.62471e-07);
    Mean->FixParameter(3,28.3609);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,1.937);
    Sigma->FixParameter(1,0.0526929);
    Sigma->FixParameter(2,-2.00484);
    Sigma->FixParameter(3,0.0204716);
    Sigma->FixParameter(4,1.17935);

  }

  if((MetaLHC18 == true) && (isDeuteron == true)){

    Mean->FixParameter(0,3.56105);
    Mean->FixParameter(1,0.000378134);
    Mean->FixParameter(2,-0.000464038);
    Mean->FixParameter(3,9.30636);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.0506275);
    Sigma->FixParameter(1,0.035813);
    Sigma->FixParameter(2,-0.0111854);
    Sigma->FixParameter(3,2.74214);
    Sigma->FixParameter(4,3);

  }

  if((MetaLHC18 == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.887204);
    Mean->FixParameter(1,0.00500025);
    Mean->FixParameter(2,-4.74017e-07);
    Mean->FixParameter(3,34.2322);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.0377498);
    Sigma->FixParameter(1,0.0586057);
    Sigma->FixParameter(2,-0.14258);
    Sigma->FixParameter(3,0.0234825);
    Sigma->FixParameter(4,34.3808);

  }

  if((MetaLHC18 == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,3.59297);
    Mean->FixParameter(1,-0.0119197);
    Mean->FixParameter(2,-0.011537);
    Mean->FixParameter(3,3.79925);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.0474376);
    Sigma->FixParameter(1,0.0362703);
    Sigma->FixParameter(2,-0.0104991);
    Sigma->FixParameter(3,2.80921);
    Sigma->FixParameter(4,3);

  }

  if((MetaLHC18 == true) && (isPion == true)){

    Mean->FixParameter(0,0.0155001);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.0055001);
    Mean->FixParameter(3,-1.06371e+09);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }

  if((MetaLHC18 == true) && (isAntiPion == true)){

    Mean->FixParameter(0,0.0181884);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.00318841);
    Mean->FixParameter(3,-3.22095e+07);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }

  // LHC18q PbPb data 
  if((LHC18q == true) && (isProton == true)){

    Mean->FixParameter(0,0.882663);
    Mean->FixParameter(1,0.00954118);
    Mean->FixParameter(2,-4.23798e-08);
    Mean->FixParameter(3,74.7086);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,-0.499897);
    Sigma->FixParameter(1,0.0498165);
    Sigma->FixParameter(2,0.423006);
    Sigma->FixParameter(3,-0.625526);
    Sigma->FixParameter(4,0.266939);

  }

  if((LHC18q == true) && (isDeuteron == true)){

    Mean->FixParameter(0,3.54);
    Mean->FixParameter(1,0.01);
    Mean->FixParameter(2,-2e-05);
    Mean->FixParameter(3,25);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,-0.0580942);
    Sigma->FixParameter(1,0.0649085);
    Sigma->FixParameter(2,0.00450883);
    Sigma->FixParameter(3,-3.21576);
    Sigma->FixParameter(4,2.42083);

  }

  if((LHC18q == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.899341);
    Mean->FixParameter(1,0.0070984);
    Mean->FixParameter(2,-0.0223373);
    Mean->FixParameter(3,1.49187);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,-0.257502);
    Sigma->FixParameter(1,0.0480432);
    Sigma->FixParameter(2,0.170963);
    Sigma->FixParameter(3,-2.64634);
    Sigma->FixParameter(4,0.265305);

  }

  if((LHC18q == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,3.54826);
    Mean->FixParameter(1,0.0177966);
    Mean->FixParameter(2,-0.00774686);
    Mean->FixParameter(3,4.20711);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.02181);
    Sigma->FixParameter(1,0.0483774);
    Sigma->FixParameter(2,-0.000824466);
    Sigma->FixParameter(3,5.53748);
    Sigma->FixParameter(4,3);

  }

  if((LHC18q == true) && (isPion == true)){

    Mean->FixParameter(0,0.0154622);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.00546224);
    Mean->FixParameter(3,-6.76792e+07);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }

  if((LHC18q == true) && (isAntiPion == true)){

    Mean->FixParameter(0,0.0155);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.0055);
    Mean->FixParameter(3,-3e+09);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }


  // LHC18r PbPb data 
  if((LHC18r == true) && (isProton == true)){

    Mean->FixParameter(0,0.886972);
    Mean->FixParameter(1,0.00847106);
    Mean->FixParameter(2,-1.94201e-08);
    Mean->FixParameter(3,93.4605);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,-0.717964);
    Sigma->FixParameter(1,0.0455849);
    Sigma->FixParameter(2,0.654092);
    Sigma->FixParameter(3,-0.213367);
    Sigma->FixParameter(4,0.378298);

  }

  if((LHC18r == true) && (isDeuteron == true)){

    Mean->FixParameter(0,3.5);
    Mean->FixParameter(1,0.03);
    Mean->FixParameter(2,-2e-05);
    Mean->FixParameter(3,25);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,-0.0424764);
    Sigma->FixParameter(1,0.0625353);
    Sigma->FixParameter(2,0.00169585);
    Sigma->FixParameter(3,-5.08514);
    Sigma->FixParameter(4,2.41159);

  }

  if((LHC18r == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.883118);
    Mean->FixParameter(1,0.00949616);
    Mean->FixParameter(2,-1.47418e-06);
    Mean->FixParameter(3,22.7463);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.103316);
    Sigma->FixParameter(1,0.0422522);
    Sigma->FixParameter(2,-0.156961);
    Sigma->FixParameter(3,0.0952412);
    Sigma->FixParameter(4,3.17989);

  }

  if((LHC18r == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,3.51);
    Mean->FixParameter(1,0.02);
    Mean->FixParameter(2,-2e-05);
    Mean->FixParameter(3,25);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,-0.0256875);
    Sigma->FixParameter(1,0.0584863);
    Sigma->FixParameter(2,0.00261934);
    Sigma->FixParameter(3,-2.12772);
    Sigma->FixParameter(4,3.38546);

  }

  if((LHC18r == true) && (isPion == true)){

    Mean->FixParameter(0,0.0154343);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.00543432);
    Mean->FixParameter(3,-3.75269e+07);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }

  if((LHC18r == true) && (isAntiPion == true)){

    Mean->FixParameter(0,0.0155);
    Mean->FixParameter(1,0);
    Mean->FixParameter(2,0.0055);
    Mean->FixParameter(3,-2.02299e+09);
    Mean->FixParameter(4,0);

    Sigma->FixParameter(0,0);
    Sigma->FixParameter(1,0.01);
    Sigma->FixParameter(2,-0.001);
    Sigma->FixParameter(3,0.1);
    Sigma->FixParameter(4,1);

  }




  double mean = Mean->Eval(pT);
  double sigma = Sigma->Eval(pT);

  Mean->Delete();
  Sigma->Delete();


  SigmaParticle = (massSq - mean)/(sigma);
  return SigmaParticle;

} // end of CalculateSigmaMassSquareTOF










// apply track cuts for protons and antiprotons
bool AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber)
{

  bool PassedParticleCuts = false;

  double Proton_pT_min, Proton_pT_max, Proton_eta_min, Proton_eta_max;
  double Proton_DCAxy_min, Proton_DCAz_min, Proton_DCAxy_max, Proton_DCAz_max;
  double Proton_TPC_RatioRowsFindableCluster_min;
  double Proton_TPC_dEdx_nSigma_max, Proton_TPC_Chi2perCluster_max, Proton_TPC_Chi2perNDF_max;
  int Proton_TPC_nCluster_min, Proton_TPC_nCrossedRows_min, Proton_TPC_nSharedCluster_max;
  double Proton_TPC_Threshold;
  double Proton_TOF_m2_nSigma_max, Proton_TOF_m2_nSigma_max_low_pTPC;
  double Proton_ITS_dEdx_nSigma_max;
  int Proton_ITS_nCluster_min;
  bool UseTOF = true;
  bool UseITS = true;


  if(fUseOpenCuts == true){

    // define open proton and antiproton track cuts
    Proton_pT_min = 0.0;
    Proton_pT_max = 4.0;
    Proton_eta_min = -0.8;
    Proton_eta_max = +0.8;
    Proton_DCAxy_min = 0.02; // cm
    Proton_DCAxy_max = 999.0; // cm
    Proton_DCAz_min = 0.02; // cm
    Proton_DCAz_max = 999.0; // cm

    Proton_TPC_RatioRowsFindableCluster_min = 0.73;
    Proton_TPC_dEdx_nSigma_max = 5.0;
    Proton_TPC_Chi2perCluster_max = 5.0;
    Proton_TPC_Chi2perNDF_max = 5.0;
    Proton_TPC_nCluster_min = 70;
    Proton_TPC_nCrossedRows_min = 60;
    Proton_TPC_nSharedCluster_max = 2;
    Proton_TPC_Threshold = 0.7;

    Proton_TOF_m2_nSigma_max = 5.0;
    Proton_TOF_m2_nSigma_max_low_pTPC = 5.0;

    Proton_ITS_dEdx_nSigma_max = 5.0;
    Proton_ITS_nCluster_min = 0;

    UseTOF = true;
    UseITS = true;

  } // end of UseOpenCuts == true



  if(fUseOpenCuts == false){

    // define closed proton and antiproton track cuts
    Proton_pT_min = 0.0;
    Proton_pT_max = 4.0;
    Proton_eta_min = -0.8;
    Proton_eta_max = +0.8;
    Proton_DCAxy_min = 0.02; // cm
    Proton_DCAxy_max = 999.0; // cm
    Proton_DCAz_min = 0.02; // cm
    Proton_DCAz_max = 999.0; // cm

    Proton_TPC_RatioRowsFindableCluster_min = 0.83;
    Proton_TPC_dEdx_nSigma_max = 3.0;
    Proton_TPC_Chi2perCluster_max = 4.0;
    Proton_TPC_Chi2perNDF_max = 4.0;
    Proton_TPC_nCluster_min = 80;
    Proton_TPC_nCrossedRows_min = 70;
    Proton_TPC_nSharedCluster_max = 0;
    Proton_TPC_Threshold = 0.7;

    Proton_TOF_m2_nSigma_max = 3.0;
    Proton_TOF_m2_nSigma_max_low_pTPC = 3.0;

    Proton_ITS_dEdx_nSigma_max = 3.0;
    Proton_ITS_nCluster_min = 1;

    UseTOF = true;
    UseITS = true;

  } // end of UseOpenCuts == false






  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;

  double p = Track.P();
  double pT = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  if(TMath::IsNaN(p)) return PassedParticleCuts;
  if(TMath::IsNaN(pT)) return PassedParticleCuts;
  if(TMath::IsNaN(pTPC)) return PassedParticleCuts;

  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Proton_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;

  // apply TPC nSigma cut
  double TPC_dEdx_nSigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
  if(TMath::IsNaN(TPC_dEdx_nSigma)) return PassedParticleCuts;
  if(TMath::Abs(TPC_dEdx_nSigma) > Proton_TPC_dEdx_nSigma_max) return PassedParticleCuts;

  // get DCA information
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  float DCAxy = xv[0];
  float DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;
  
  // apply DCAxy cut
  if((TMath::Abs(DCAxy) < Proton_DCAxy_min) || (TMath::Abs(DCAxy) > Proton_DCAxy_max)) return PassedParticleCuts;

  // apply DCAz cut
  if((TMath::Abs(DCAz) < Proton_DCAz_min) || (TMath::Abs(DCAz) > Proton_DCAz_max)) return PassedParticleCuts;

  // apply pT cut
  if(pT < Proton_pT_min || pT > Proton_pT_max) return PassedParticleCuts;

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(eta < Proton_eta_min || eta > Proton_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  int TPC_nCluster = Track.GetNcls(1);
  if(TPC_nCluster < Proton_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  int TPC_nCrossedRows = Track.GetTPCCrossedRows();
  if(TPC_nCrossedRows < Proton_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply shared cluster cut for TPC
  int TPC_nSharedCluster = Track.GetTPCnclsS();
  if(TPC_nSharedCluster > Proton_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  int TPC_nFindableCluster = Track.GetTPCNclsF();
  double TPC_RatioRowsFindableCluster = -999.0;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((double)TPC_nCrossedRows / (double)TPC_nFindableCluster);
  if(TPC_RatioRowsFindableCluster < Proton_TPC_RatioRowsFindableCluster_min) return PassedParticleCuts;

  // receive TPC chi2
  double TPC_Chi2 = Track.GetTPCchi2();
  if(TMath::IsNaN(TPC_Chi2)) return PassedParticleCuts;

  // compute TPC ndf
  double TPC_NDF = -999.0;
  if(TPC_nCluster > 5) TPC_NDF = (2*TPC_nCluster-5);

  // apply TPC chi2 per cluster cut
  double TPC_Chi2perCluster = TPC_Chi2/TPC_nCluster;
  if(TPC_Chi2perCluster > Proton_TPC_Chi2perCluster_max) return PassedParticleCuts;

  // apply TPC chi2 per ndf cut
  double TPC_Chi2perNDF = TPC_Chi2/TPC_NDF;
  if(TPC_Chi2perNDF > Proton_TPC_Chi2perNDF_max) return PassedParticleCuts; 




  // check if ITS information is available
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  bool ITSisOK = false;
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;


  if((ITSisOK == true) && (isMatter == true)) h_Proton_ITS_dEdx_NoTOFcutNoITScut->Fill(p,Track.GetITSsignal());
  if((ITSisOK == true) && (isMatter == false)) h_AntiProton_ITS_dEdx_NoTOFcutNoITScut->Fill(p,Track.GetITSsignal());


  if((ITSisOK == true) && (UseITS == true)){

    int ParticleSpecies = 0;
    if(isMatter) ParticleSpecies = 1;
    if(!isMatter) ParticleSpecies = 3;
    
    double ITS_dEdx_Sigma = CalculateSigmadEdxITS(Track,ParticleSpecies,RunNumber);
    if(TMath::Abs(ITS_dEdx_Sigma) > Proton_ITS_dEdx_nSigma_max) return PassedParticleCuts;

    // apply ITS cluster cut
    double nClusterITS = Track.GetITSNcls();
    if(TMath::IsNaN(nClusterITS)) return PassedParticleCuts;
    if(nClusterITS < Proton_ITS_nCluster_min) return PassedParticleCuts;

  } // end of ITSisOK




  if((TOFisOK == true) && (isMatter == true)) h_Proton_TOF_m2_NoTOFcut->Fill(pT,CalculateMassSquareTOF(Track));
  if((TOFisOK == true) && (isMatter == false)) h_AntiProton_TOF_m2_NoTOFcut->Fill(pT,CalculateMassSquareTOF(Track));

  if((TOFisOK == true) && (UseTOF == true)){

    int ParticleSpecies = 0;
    if(isMatter) ParticleSpecies = 1;
    if(!isMatter) ParticleSpecies = 3;

    double TOF_m2	  = CalculateMassSquareTOF(Track);
    double TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(pT,TOF_m2,ParticleSpecies,RunNumber);

    // apply tight TOF m2 cut above pTPC threshold   
    if(pTPC >= Proton_TPC_Threshold){

      if(TMath::Abs(TOF_m2_nSigma) > Proton_TOF_m2_nSigma_max) return PassedParticleCuts;

    }

    // apply loose TOF m2 cut below pTPC threshold
    if(pTPC < Proton_TPC_Threshold){

      if(TMath::Abs(TOF_m2_nSigma) > Proton_TOF_m2_nSigma_max_low_pTPC) return PassedParticleCuts;

    }
 
  } // end of TOFisOK





  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckProtonCuts







// apply track cuts for Pions and antiPions
bool AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckPionCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, double Pion_DCAxy_min, double Pion_DCAz_min)
{

  bool PassedParticleCuts = false;

  double Pion_pT_min, Pion_pT_max, Pion_eta_min, Pion_eta_max;
  double Pion_DCAxy_max, Pion_DCAz_max;
  double Pion_TPC_RatioRowsFindableCluster_min;
  double Pion_TPC_dEdx_nSigma_max, Pion_TPC_Chi2perCluster_max, Pion_TPC_Chi2perNDF_max;
  int Pion_TPC_nCluster_min, Pion_TPC_nCrossedRows_min, Pion_TPC_nSharedCluster_max;
  double Pion_TPC_Threshold;
  double Pion_TOF_m2_nSigma_max, Pion_TOF_m2_nSigma_max_low_pTPC;
  double Pion_ITS_dEdx_nSigma_max;
  int Pion_ITS_nCluster_min;
  bool UseTOF = false;
  bool UseITS = false;

  if(fUseOpenCuts == true){

    // define open Pion and antiPion track cuts
    Pion_pT_min = 0.0;
    Pion_pT_max = 4.0;
    Pion_eta_min = -0.8;
    Pion_eta_max = +0.8;
    Pion_DCAxy_max = 999.0; // cm
    Pion_DCAz_max = 999.0; // cm

    Pion_TPC_RatioRowsFindableCluster_min = 0.73;
    Pion_TPC_dEdx_nSigma_max = 5.0;
    Pion_TPC_Chi2perCluster_max = 5.0;
    Pion_TPC_Chi2perNDF_max = 5.0;
    Pion_TPC_nCluster_min = 60;
    Pion_TPC_nCrossedRows_min = 50;
    Pion_TPC_nSharedCluster_max = 2;
    Pion_TPC_Threshold = 4.0;

    Pion_TOF_m2_nSigma_max = 5.0;
    Pion_TOF_m2_nSigma_max_low_pTPC = 5.0;

    Pion_ITS_dEdx_nSigma_max = 5.0;
    Pion_ITS_nCluster_min = 0;

    UseTOF = true;
    UseITS = true;

  } // end of UseOpenCuts == true



  if(fUseOpenCuts == false){

    // define closed Pion and antiPion track cuts
    Pion_pT_min = 0.0;
    Pion_pT_max = 4.0;
    Pion_eta_min = -0.8;
    Pion_eta_max = +0.8;
    Pion_DCAxy_max = 999.0; // cm
    Pion_DCAz_max = 999.0; // cm

    Pion_TPC_RatioRowsFindableCluster_min = 0.83;
    Pion_TPC_dEdx_nSigma_max = 3.0;
    Pion_TPC_Chi2perCluster_max = 4.0;
    Pion_TPC_Chi2perNDF_max = 4.0;
    Pion_TPC_nCluster_min = 80;
    Pion_TPC_nCrossedRows_min = 70;
    Pion_TPC_nSharedCluster_max = 0;
    Pion_TPC_Threshold = 0.7;

    Pion_TOF_m2_nSigma_max = 3.0;
    Pion_TOF_m2_nSigma_max_low_pTPC = 3.0;

    Pion_ITS_dEdx_nSigma_max = 3.0;
    Pion_ITS_nCluster_min = 2;

    UseTOF = true;
    UseITS = true;

  } // end of UseOpenCuts == false



  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;

  double p = Track.P();
  double pT = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  if(TMath::IsNaN(p)) return PassedParticleCuts;
  if(TMath::IsNaN(pT)) return PassedParticleCuts;
  if(TMath::IsNaN(pTPC)) return PassedParticleCuts;

  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Pion_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;

  // apply TPC nSigma cut
  double TPC_dEdx_nSigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kPion);
  if(TMath::IsNaN(TPC_dEdx_nSigma)) return PassedParticleCuts;
  if(TMath::Abs(TPC_dEdx_nSigma) > Pion_TPC_dEdx_nSigma_max) return PassedParticleCuts;

  // get DCA information
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  float DCAxy = xv[0];
  float DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;
  
  // apply DCAxy cut
  if((TMath::Abs(DCAxy) < Pion_DCAxy_min) || (TMath::Abs(DCAxy) > Pion_DCAxy_max)) return PassedParticleCuts;

  // apply DCAz cut
  if((TMath::Abs(DCAz) < Pion_DCAz_min) || (TMath::Abs(DCAz) > Pion_DCAz_max)) return PassedParticleCuts;

  // apply pT cut
  if(pT < Pion_pT_min || pT > Pion_pT_max) return PassedParticleCuts;

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(eta < Pion_eta_min || eta > Pion_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  int TPC_nCluster = Track.GetNcls(1);
  if(TPC_nCluster < Pion_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  int TPC_nCrossedRows = Track.GetTPCCrossedRows();
  if(TPC_nCrossedRows < Pion_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply shared cluster cut for TPC
  int TPC_nSharedCluster = Track.GetTPCnclsS();
  if(TPC_nSharedCluster > Pion_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  int TPC_nFindableCluster = Track.GetTPCNclsF();
  double TPC_RatioRowsFindableCluster = -999.0;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((double)TPC_nCrossedRows / (double)TPC_nFindableCluster);
  if(TPC_RatioRowsFindableCluster < Pion_TPC_RatioRowsFindableCluster_min) return PassedParticleCuts;

  // receive TPC chi2
  double TPC_Chi2 = Track.GetTPCchi2();
  if(TMath::IsNaN(TPC_Chi2)) return PassedParticleCuts;

  // compute TPC ndf
  double TPC_NDF = -999.0;
  if(TPC_nCluster > 5) TPC_NDF = (2*TPC_nCluster-5);

  // apply TPC chi2 per cluster cut
  double TPC_Chi2perCluster = TPC_Chi2/TPC_nCluster;
  if(TPC_Chi2perCluster > Pion_TPC_Chi2perCluster_max) return PassedParticleCuts;

  // apply TPC chi2 per ndf cut
  double TPC_Chi2perNDF = TPC_Chi2/TPC_NDF;
  if(TPC_Chi2perNDF > Pion_TPC_Chi2perNDF_max) return PassedParticleCuts; 




  // check if ITS information is available
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  bool ITSisOK = false;
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;


  if((ITSisOK == true) && (isMatter == true))	h_Pion_ITS_dEdx_NoTOFcutNoITScut->Fill(p,Track.GetITSsignal());
  if((ITSisOK == true) && (isMatter == false))	h_AntiPion_ITS_dEdx_NoTOFcutNoITScut->Fill(p,Track.GetITSsignal());


  if((ITSisOK == true) && (UseITS == true)){

    int ParticleSpecies = 0;
    if(isMatter) ParticleSpecies = 5;
    if(!isMatter) ParticleSpecies = 6;
  
    double ITS_dEdx_Sigma = CalculateSigmadEdxITS(Track,ParticleSpecies,RunNumber);
    if(TMath::Abs(ITS_dEdx_Sigma) > Pion_ITS_dEdx_nSigma_max) return PassedParticleCuts;

    // apply ITS cluster cut
    double nClusterITS = Track.GetITSNcls();
    if(TMath::IsNaN(nClusterITS)) return PassedParticleCuts;
    if(nClusterITS < Pion_ITS_nCluster_min) return PassedParticleCuts;

  } // end of ITSisOK




  if((TOFisOK == true) && (isMatter == true))	h_Pion_TOF_m2_NoTOFcut->Fill(pT,CalculateMassSquareTOF(Track));
  if((TOFisOK == true) && (isMatter == false))	h_AntiPion_TOF_m2_NoTOFcut->Fill(pT,CalculateMassSquareTOF(Track));

  if((TOFisOK == true) && (UseTOF == true)){

    int ParticleSpecies = 0;
    if(isMatter) ParticleSpecies = 5;
    if(!isMatter) ParticleSpecies = 6;

    double TOF_m2	  = CalculateMassSquareTOF(Track);
    double TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(pT,TOF_m2,ParticleSpecies,RunNumber);

    // apply tight TOF m2 cut above pTPC threshold   
    if(pTPC >= Pion_TPC_Threshold){

      if(TMath::Abs(TOF_m2_nSigma) > Pion_TOF_m2_nSigma_max) return PassedParticleCuts;

    }

    // apply loose TOF m2 cut below pTPC threshold
    if(pTPC < Pion_TPC_Threshold){

      if(TMath::Abs(TOF_m2_nSigma) > Pion_TOF_m2_nSigma_max_low_pTPC) return PassedParticleCuts;

    }

  } // end of TOFisOK





  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckPionCuts






































bool AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckLambdaCuts(AliAODv0 &v0, double PrimaryVertexPos[3], AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, double MassInvariant, double WrongMassInvariant)
{

  bool PassedParticleCuts = false;


  double eta_max = 0.8;
  double pT_min = 0.0; // GeV/c
  double pT_max = 999.0; // GeV/c
  double DecayRadius_min = 0.1; // cm
  double DecayRadius_max = 100.0; // cm
  double LambdaMassVariation_max = 0.01; // GeV/c²
  double DCAv0ToPrimaryVertex_max = 1.5; // cm
  double CosinePointingAngle_min = 0.997;

  bool UseReconstructionOnTheFly = true;

  double DCAv0Daughters_min = 0.0;
  double DCAv0Daughters_max = 2.0;
 
 
  bool IsReconstructedOnTheFly = v0.GetOnFlyStatus();
  if(!(IsReconstructedOnTheFly == UseReconstructionOnTheFly)) return PassedParticleCuts;
 
  int nDaughters = v0.GetNDaughters();
  if(!(nDaughters == 2)) return PassedParticleCuts;

  int nProngs = v0.GetNProngs();
  if(!(nProngs == 2)) return PassedParticleCuts;

  int Charge = v0.GetCharge();
  if(!(Charge == 0)) return PassedParticleCuts;

  int Charge1 = v0.ChargeProng(0);
  int Charge2 = v0.ChargeProng(1);
  if(Charge1 == Charge2) return PassedParticleCuts;
 
  double eta = v0.PseudoRapV0();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(TMath::Abs(eta) > eta_max) return PassedParticleCuts;

  double DCAv0ToPrimaryVertex = v0.DcaV0ToPrimVertex();
  if(TMath::IsNaN(DCAv0ToPrimaryVertex)) return PassedParticleCuts;
  if(TMath::Abs(DCAv0ToPrimaryVertex) > DCAv0ToPrimaryVertex_max) return PassedParticleCuts;

  double DecayRadius = v0.RadiusV0();
  if(TMath::IsNaN(DecayRadius)) return PassedParticleCuts;
  if((TMath::Abs(DecayRadius) < DecayRadius_min) || (TMath::Abs(DecayRadius) > DecayRadius_max)) return PassedParticleCuts;

  double px = v0.MomV0X();
  double py = v0.MomV0Y();
  if(TMath::IsNaN(px)) return PassedParticleCuts;
  if(TMath::IsNaN(py)) return PassedParticleCuts;
  double pT = TMath::Sqrt((px * px) + (py * py));
  if((pT < pT_min) || (pT > pT_max)) return PassedParticleCuts; 

  double CosinePointingAngle = v0.CosPointingAngle(PrimaryVertexPos);
  if(TMath::IsNaN(CosinePointingAngle)) return PassedParticleCuts;
  if(TMath::Abs(CosinePointingAngle) < CosinePointingAngle_min) return PassedParticleCuts; 

  double DCAv0Daughters = v0.DcaV0Daughters();
  if(TMath::IsNaN(DCAv0Daughters)) return PassedParticleCuts;
  if((TMath::Abs(DCAv0Daughters) < DCAv0Daughters_min) || (TMath::Abs(DCAv0Daughters) > DCAv0Daughters_max)) return PassedParticleCuts;




  if(fUseOpenCuts == true){

    if(MassInvariant < 1.07 || MassInvariant > 1.17) return PassedParticleCuts;
    if(WrongMassInvariant < 1.4) return PassedParticleCuts;

  }

  if(fUseOpenCuts == false){

    const double LambdaMassPDG = 1.115683; // GeV/c²
    if(MassInvariant < (LambdaMassPDG-LambdaMassVariation_max) || MassInvariant > (LambdaMassPDG+LambdaMassVariation_max)) return PassedParticleCuts;
    if(WrongMassInvariant < 1.4) return PassedParticleCuts;

  }
 

  PassedParticleCuts = true;
  return PassedParticleCuts;

}







// apply track cuts for deuterons and antideuterons
bool AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber)
{

  bool PassedParticleCuts = false;


  // define open deuteron and antideuteron track cuts
  double Deuteron_pT_min = 0.5;
  double Deuteron_pT_max = 3.0;
  double Deuteron_eta_min = -0.8;
  double Deuteron_eta_max = +0.8;
  double Deuteron_DCAxy_max = 0.3; // cm
  double Deuteron_DCAz_max = 0.2; // cm

  double Deuteron_TPC_RatioRowsFindableCluster_min = 0.73;
  double Deuteron_TPC_dEdx_nSigma_max = 4.0;
  double Deuteron_TPC_Chi2perCluster_max = 5.0;
  double Deuteron_TPC_Chi2perNDF_max = 5.0;
  int Deuteron_TPC_nCluster_min = 70;
  int Deuteron_TPC_nCrossedRows_min = 60;
  int Deuteron_TPC_nSharedCluster_max = 2;
  double Deuteron_TPC_Threshold = 1.2;

  double Deuteron_TOF_m2_nSigma_max = 4.0;

  double Deuteron_ITS_dEdx_nSigma_max = 4.0;
  int Deuteron_ITS_nCluster_min = 1;

  bool UseTOF = true;
  bool UseITS = true;
  bool Extend_pT_range = true;

  double Pion_TPC_dEdx_nSigma_max     = 3.0;
  double Kaon_TPC_dEdx_nSigma_max     = 3.0;
  double Proton_TPC_dEdx_nSigma_max   = 3.0;
  double Electron_TPC_dEdx_nSigma_max = 3.0;
  double Muon_TPC_dEdx_nSigma_max     = 3.0;






  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;

  double p = Track.P();
  double pT = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  if(TMath::IsNaN(p)) return PassedParticleCuts;
  if(TMath::IsNaN(pT)) return PassedParticleCuts;
  if(TMath::IsNaN(pTPC)) return PassedParticleCuts;

  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Deuteron_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;


  // apply TPC nSigma cut
  double TPC_dEdx_nSigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
  if(TMath::IsNaN(TPC_dEdx_nSigma)) return PassedParticleCuts;

  if(fUseOpenCuts == false){
    if(TMath::Abs(TPC_dEdx_nSigma) > Deuteron_TPC_dEdx_nSigma_max) return PassedParticleCuts;
  }

  if(fUseOpenCuts == true){
    if(TMath::Abs(TPC_dEdx_nSigma) > 8.0) return PassedParticleCuts;
  }

  // get DCA information
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  float DCAxy = xv[0];
  float DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;
  
  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Deuteron_DCAxy_max) return PassedParticleCuts;

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Deuteron_DCAz_max) return PassedParticleCuts;

  // apply pT cut
  if(pT < Deuteron_pT_min || pT > Deuteron_pT_max) return PassedParticleCuts;

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(eta < Deuteron_eta_min || eta > Deuteron_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  int TPC_nCluster = Track.GetNcls(1);
  if(TPC_nCluster < Deuteron_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  int TPC_nCrossedRows = Track.GetTPCCrossedRows();
  if(TPC_nCrossedRows < Deuteron_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply zero shared cluster cut for TPC
  int TPC_nSharedCluster = Track.GetTPCnclsS();
  if(TPC_nSharedCluster > Deuteron_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  int TPC_nFindableCluster = Track.GetTPCNclsF();
  double TPC_RatioRowsFindableCluster = -999.0;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((double)TPC_nCrossedRows / (double)TPC_nFindableCluster);
  if(TPC_RatioRowsFindableCluster < Deuteron_TPC_RatioRowsFindableCluster_min) return PassedParticleCuts;

  // receive TPC chi2
  double TPC_Chi2 = Track.GetTPCchi2();
  if(TMath::IsNaN(TPC_Chi2)) return PassedParticleCuts;

  // compute TPC ndf
  double TPC_NDF = -999.0;
  if(TPC_nCluster > 5) TPC_NDF = (2*TPC_nCluster-5);

  // apply TPC chi2 per cluster cut
  double TPC_Chi2perCluster = TPC_Chi2/TPC_nCluster;
  if(TPC_Chi2perCluster > Deuteron_TPC_Chi2perCluster_max) return PassedParticleCuts;

  // apply TPC chi2 per ndf cut
  double TPC_Chi2perNDF = TPC_Chi2/TPC_NDF;
  if(TPC_Chi2perNDF > Deuteron_TPC_Chi2perNDF_max) return PassedParticleCuts; 


  if(Extend_pT_range == true){

    // cut out dEdx band of other particles above pTPC = 1.6 GeV/c²
    double TPC_dEdx_nSigma_Pion = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kPion);
    if(TMath::IsNaN(TPC_dEdx_nSigma_Pion)) return PassedParticleCuts;
    if((pT >= 1.6) && (TMath::Abs(TPC_dEdx_nSigma_Pion) < Pion_TPC_dEdx_nSigma_max)) return PassedParticleCuts;

    double TPC_dEdx_nSigma_Kaon = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kKaon);
    if(TMath::IsNaN(TPC_dEdx_nSigma_Kaon)) return PassedParticleCuts;
    if((pT >= 1.6) && (TMath::Abs(TPC_dEdx_nSigma_Kaon) < Kaon_TPC_dEdx_nSigma_max)) return PassedParticleCuts;

    double TPC_dEdx_nSigma_Proton = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
    if(TMath::IsNaN(TPC_dEdx_nSigma_Proton)) return PassedParticleCuts;
    if((pT >= 1.6) && (TMath::Abs(TPC_dEdx_nSigma_Proton) < Proton_TPC_dEdx_nSigma_max)) return PassedParticleCuts;

    double TPC_dEdx_nSigma_Electron = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kElectron);
    if(TMath::IsNaN(TPC_dEdx_nSigma_Electron)) return PassedParticleCuts;
    if((pT >= 1.6) && (TMath::Abs(TPC_dEdx_nSigma_Electron) < Electron_TPC_dEdx_nSigma_max)) return PassedParticleCuts;

    double TPC_dEdx_nSigma_Muon = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kMuon);
    if(TMath::IsNaN(TPC_dEdx_nSigma_Muon)) return PassedParticleCuts;
    if((pT >= 1.6) && (TMath::Abs(TPC_dEdx_nSigma_Muon) < Muon_TPC_dEdx_nSigma_max)) return PassedParticleCuts;

  } // end of Extend_pT_range


  // check if ITS information is available
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  bool ITSisOK = false;
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;


  if((ITSisOK == true) && (isMatter == true)) h_Deuteron_ITS_dEdx_NoTOFcutNoITScut->Fill(p,Track.GetITSsignal());
  if((ITSisOK == true) && (isMatter == false)) h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut->Fill(p,Track.GetITSsignal());


  if((ITSisOK == true) && (UseITS == true)){

    int ParticleSpecies = 0;
    if(isMatter) ParticleSpecies = 2;
    if(!isMatter) ParticleSpecies = 4;
    
    double ITS_dEdx_Sigma = CalculateSigmadEdxITS(Track,ParticleSpecies,RunNumber);
    if(TMath::Abs(ITS_dEdx_Sigma) > Deuteron_ITS_dEdx_nSigma_max) return PassedParticleCuts;

    // apply ITS cluster cut
    double nClusterITS = Track.GetITSNcls();
    if(TMath::IsNaN(nClusterITS)) return PassedParticleCuts;
    if(nClusterITS < Deuteron_ITS_nCluster_min) return PassedParticleCuts;

  } // end of ITSisOK





  if((TOFisOK == true) && (isMatter == true)) h_Deuteron_TOF_m2_NoTOFcut->Fill(pT,CalculateMassSquareTOF(Track));
  if((TOFisOK == true) && (isMatter == false)) h_AntiDeuteron_TOF_m2_NoTOFcut->Fill(pT,CalculateMassSquareTOF(Track));


  if((TOFisOK == true) && (UseTOF == true)){

    int ParticleSpecies = 0;
    if(isMatter) ParticleSpecies = 2;
    if(!isMatter) ParticleSpecies = 4;

    double TOF_m2	  = CalculateMassSquareTOF(Track);
    double TOF_m2_nSigma  = CalculateSigmaMassSquareTOF(pT,TOF_m2,ParticleSpecies,RunNumber);

    if(TMath::Abs(TOF_m2_nSigma) > Deuteron_TOF_m2_nSigma_max) return PassedParticleCuts;

  } // end of TOFisOK





  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckDeuteronCuts
















double AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateSigmadEdxITS(AliAODTrack &Track, int ParticleSpecies, int RunNumber){

  double SigmaParticle = -999.0;
  double SignalITS = Track.GetITSsignal();
  if(TMath::IsNaN(SignalITS)) return SigmaParticle;

  bool MetaLHC16 = false;
  bool MetaLHC17 = false;
  bool MetaLHC18 = false;
  bool LHC18q = false;
  bool LHC18r = false;
  bool LHC20g7a = false;
  bool LHC20g7b = false;
  bool LHC22f3 = false;

  if(fIsMC == false){
  
    if((RunNumber >= 252235) && (RunNumber <= 264347)) MetaLHC16 = true;
    if((RunNumber >= 270581) && (RunNumber <= 282704)) MetaLHC17 = true;
    if((RunNumber >= 285009) && (RunNumber <= 294925)) MetaLHC18 = true;
    if((RunNumber >= 295585) && (RunNumber <= 296623)) LHC18q = true;
    if((RunNumber >= 296690) && (RunNumber <= 297585)) LHC18r = true;

  } // end of fIsMC == false


  if(fIsMC == true){

    if(fCollisionSystem == 1) LHC20g7a	= true;
    if(fCollisionSystem == 2) LHC20g7b	= true;
    if(fCollisionSystem > 2)  LHC22f3	= true;

  } // end of fIsMC == true


  bool isProton	      = false;
  bool isDeuteron     = false;
  bool isAntiProton   = false;
  bool isAntiDeuteron = false;
  bool isPion	      = false;
  bool isAntiPion     = false;

  if(ParticleSpecies == 1) isProton = true;
  if(ParticleSpecies == 2) isDeuteron = true;
  if(ParticleSpecies == 3) isAntiProton = true;
  if(ParticleSpecies == 4) isAntiDeuteron = true;
  if(ParticleSpecies == 5) isPion = true;
  if(ParticleSpecies == 6) isAntiPion = true;

  double p = Track.P();
  double Mass = 0.0;
  if((isProton == true) || (isAntiProton == true))	Mass = AliPID::ParticleMass(AliPID::kProton);
  if((isDeuteron == true) || (isAntiDeuteron == true))	Mass = AliPID::ParticleMass(AliPID::kDeuteron);
  if((isPion == true) || (isAntiPion == true))		Mass = AliPID::ParticleMass(AliPID::kPion);


  TF1 *Mean = new TF1("Mean","[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.0,6.0);
  Mean->FixParameter(5,1);
  Mean->FixParameter(6,Mass);


  // LHC20g7a
  if((LHC20g7a == true) && (isProton == true)){
    
    Mean->FixParameter(0,5.25076e-18);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,127.35);
    Mean->FixParameter(4,9388.68);

  }

  if((LHC20g7a == true) && (isDeuteron == true)){

    Mean->FixParameter(0,5.13883e-18);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,9703.93);

  }

  if(LHC20g7a == true && isAntiProton == true){
    
    Mean->FixParameter(0,5.45646e-27);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,127.35);
    Mean->FixParameter(4,6792.38);

  }

  if(LHC20g7a == true && isAntiDeuteron == true){

    Mean->FixParameter(0,6.65152e-21);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,8811.37);

  }

  if(LHC20g7a == true && isPion == true){

    Mean->FixParameter(0,1.8487e-09);
    Mean->FixParameter(1,-163771);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,15658.8);

  }

  if(LHC20g7a == true && isAntiPion == true){

    Mean->FixParameter(0,3.14614e-10);
    Mean->FixParameter(1,-157437);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,14760.3);

  }





  // LHC20g7b
  if((LHC20g7b == true) && (isProton == true)){
    
    Mean->FixParameter(0,1.29354e-13);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,127.757);
    Mean->FixParameter(4,11475.6);
    
  }

  if((LHC20g7b == true) && (isDeuteron == true)){

    Mean->FixParameter(0,3.34314e-18);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,9503.91);

  }

  if((LHC20g7b == true) && (isAntiProton == true)){
    
    Mean->FixParameter(0,4.10675e-08);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,9.48982e+12);
    Mean->FixParameter(4,16775.7);

  }

  if((LHC20g7b == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,3.34314e-18);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,9503.91);

  }
 
  if((LHC20g7b == true) && (isPion == true)){

    Mean->FixParameter(0,1.16208e-06);
    Mean->FixParameter(1,-121114);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,19595);

  }

  if((LHC20g7b == true) && (isAntiPion == true)){

    Mean->FixParameter(0,2.0497);
    Mean->FixParameter(1,-123189);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,47257.2);

  }




  // copied from data
  if((LHC22f3 == true) && ((isProton == true) || (isAntiProton == true))){

    Mean->FixParameter(0,2.36861e-07);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,9.55834);
    Mean->FixParameter(4,17081);

  }

  if((LHC22f3 == true) && (isDeuteron == true)){

    Mean->FixParameter(0,8.35954e-20);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,8627.27);
    Mean->FixParameter(4,8861.17);

  }


  if((LHC22f3 == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,8.57849e-21);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,8488.53);

  }

  if((LHC22f3 == true) && (isPion == true)){

    Mean->FixParameter(0,1.08232e-08);
    Mean->FixParameter(1,-84115.4);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,17225.4);

  }

  if((LHC22f3 == true) && (isAntiPion == true)){

    Mean->FixParameter(0,1.20941e-08);
    Mean->FixParameter(1,-84115.4);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,899.867);
    Mean->FixParameter(4,16967.4);

  }



  if((fIsMC == false) && ((isProton == true) || (isAntiProton == true))){

    Mean->FixParameter(0,2.36861e-07);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,9.55834);
    Mean->FixParameter(4,17081);

  }

  if((fIsMC == false) && ((isDeuteron == true) || (isAntiDeuteron == true))){

    Mean->FixParameter(0,7.41722e-06);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,11249.3);
    Mean->FixParameter(4,19828.9);

  }

  if((fIsMC == false) && ((isPion == true) || (isAntiPion == true))){

    Mean->FixParameter(0,2.02983e-12);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,1187.83);
    Mean->FixParameter(4,12309.8);

  }



  double mean = Mean->Eval(p);
  Mean->Delete();

  double Resolution = 0.0;

  if(((isProton == true) || (isAntiProton == true))	&& (MetaLHC16 == true)) Resolution = 0.126;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && (MetaLHC16 == true)) Resolution = 9.14588e-02;
  if(((isPion == true) || (isAntiPion == true)) && (MetaLHC16 == true)) Resolution = 1.28499e-01;

  if(((isProton == true) || (isAntiProton == true))	&& (MetaLHC17 == true)) Resolution = 1.34216e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && (MetaLHC17 == true)) Resolution = 9.00246e-02;
  if(((isPion == true) || (isAntiPion == true)) && (MetaLHC17 == true)) Resolution = 1.31156e-01;

  if(((isProton == true) || (isAntiProton == true))	&& (MetaLHC18 == true)) Resolution = 1.32506e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && (MetaLHC18 == true)) Resolution = 8.82121e-02;
  if(((isPion == true) || (isAntiPion == true)) && (MetaLHC18 == true)) Resolution = 1.32900e-01;

  if(((isProton == true) || (isAntiProton == true))	&& ((LHC18q == true) || (LHC18r == true))) Resolution = 0.10;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && ((LHC18q == true) || (LHC18r == true))) Resolution = 0.10;
  if(((isPion == true) || (isAntiPion == true)) && ((LHC18q == true) || (LHC18r == true))) Resolution = 1.71541e-01;

  if(((isProton == true) || (isAntiProton == true))	&& (LHC20g7a == true)) Resolution = 1.31668e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true))	&& (LHC20g7a == true)) Resolution = 9.46937e-02;
  if(((isPion == true) || (isAntiPion == true))	&& (LHC20g7a == true))	       Resolution = 1.73629e-01;

  if(((isProton == true) || (isAntiProton == true))	&& (LHC20g7b == true)) Resolution = 1.30878e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true))	&& (LHC20g7b == true)) Resolution = 9.46815e-02;
  if(((isPion == true) || (isAntiPion == true))	&& (LHC20g7b == true)) Resolution = 1.06914e-01;

  if(((isProton == true) || (isAntiProton == true))	&& (LHC22f3 == true)) Resolution = 1.10359e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true))	&& (LHC22f3 == true)) Resolution = 9.35349e-02;
  if(((isPion == true) || (isAntiPion == true))	&& (LHC22f3 == true)) Resolution = 8.22958e-02;

  double ScaleFactor = 1.0-(Resolution);
  double sigma = (mean*ScaleFactor) - mean;
  if(TMath::Abs(sigma) < 0.0001) return -999.0;

  SigmaParticle = (mean - SignalITS) / (sigma);

  return SigmaParticle;

} // end of CalculateSigmadEdxITS






double AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateInvariantMassLambda(double Momentum1[3], double Momentum2[3], bool isMatter){


  const double MassProton   = 0.938272088;    // GeV/c²
  const double MassPion     = 0.13957039;     // GeV/c²

  double m1, m2;
  double p1x, p1y, p1z;
  double p2x, p2y, p2z;

  bool CalculateLambda = false;
  bool CalculateAntiLambda = false;

  if(isMatter == true)	CalculateLambda = true;
  if(isMatter == false)	CalculateAntiLambda = true;

  if(CalculateLambda == true){

    m1 = MassProton;
    p1x = Momentum1[0];
    p1y = Momentum1[1];
    p1z = Momentum1[2];

    m2 = MassPion;
    p2x = Momentum2[0];
    p2y = Momentum2[1];
    p2z = Momentum2[2];

  }

  if(CalculateAntiLambda == true){

    m1 = MassPion;
    p1x = Momentum1[0];
    p1y = Momentum1[1];
    p1z = Momentum1[2];

    m2 = MassProton;
    p2x = Momentum2[0];
    p2y = Momentum2[1];
    p2z = Momentum2[2];

  }

  double E1 = TMath::Sqrt((m1*m1) + (p1x*p1x) + (p1y*p1y) + (p1z*p1z));
  double E2 = TMath::Sqrt((m2*m2) + (p2x*p2x) + (p2y*p2y) + (p2z*p2z)); 

  double MassInvariant = TMath::Sqrt(((E1+E2)*(E1+E2)) - ((p1x+p2x)*(p1x+p2x)) - ((p1y+p2y)*(p1y+p2y)) - ((p1z+p2z)*(p1z+p2z)));

  return MassInvariant;

} // end of CalculateInvariantMassLambda










