#include "TChain.h" 
#include "TTree.h"
#include "TList.h"
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
  fLambda_px(-999.0),
  fLambda_py(-999.0),
  fLambda_pz(-999.0),
  fLambda_px_Generated(-999.0),
  fLambda_py_Generated(-999.0),
  fLambda_pz_Generated(-999.0),
  fLambda_Eta(-999.0),
  fLambda_Phi(-999.0),
  fLambda_TransverseRadius(-999.0),
  fLambda_CosinePointingAngle(-999.0),
  fLambda_DCAv0ToPrimaryVertex(-999.0),
  fLambda_DCAv0Daughters(-999.0),
  fLambda_Alpha(-999.0),
  fLambda_qT(-999.0),
  fLambda_DecayLength(-999.0),
  fLambda_PDG_Daughter1(0),
  fLambda_PDG_Daughter2(0),
  fLambda_PDG_v01(0),
  fLambda_PDG_v02(0),
  fLambda_PDG_Mother1(0),
  fLambda_PDG_Mother2(0),
  fLambda_SameV0(0),
  fLambda_Event_Centrality(-999.0),
  fLambda_Event_PrimaryVertexZ(-999.0),
  fLambda_Event_BField(0),
  fLambda_Event_Multiplicity(0),
  fLambda_Event_Identifier(0),
  fLambda_Event_RunNumber(0),
  fLambda_Event_IsFirstParticle(0),
  fLambda_Daughter_Proton_px(-999.0),
  fLambda_Daughter_Proton_py(-999.0),
  fLambda_Daughter_Proton_pz(-999.0),
  fLambda_Daughter_Proton_px_Generated(-999.0),
  fLambda_Daughter_Proton_py_Generated(-999.0),
  fLambda_Daughter_Proton_pz_Generated(-999.0),
  fLambda_Daughter_Proton_px_DecayVertex(-999.0),
  fLambda_Daughter_Proton_py_DecayVertex(-999.0),
  fLambda_Daughter_Proton_pz_DecayVertex(-999.0),
  fLambda_Daughter_Proton_pTPC(-999.0),
  fLambda_Daughter_Proton_Eta(-999.0),
  fLambda_Daughter_Proton_Phi(-999.0),
  fLambda_Daughter_Proton_TPC_Chi2(-999.0),
  fLambda_Daughter_Proton_TPC_dEdx(-999.0),
  fLambda_Daughter_Proton_TPC_dEdx_Sigma(-999.0),
  fLambda_Daughter_Proton_TOF_Mass2(-999.0),
  fLambda_Daughter_Proton_TOF_Mass2_Sigma(-999.0),
  fLambda_Daughter_Proton_ITS_dEdx(-999.0),
  fLambda_Daughter_Proton_ITS_dEdx_Sigma(-999.0),
  fLambda_Daughter_Proton_DCAxy(-999.0),
  fLambda_Daughter_Proton_DCAz(-999.0),
  fLambda_Daughter_Proton_TPC_nCrossedRows(0),
  fLambda_Daughter_Proton_TPC_nFindableCluster(0),
  fLambda_Daughter_Proton_TPC_nCluster(0),
  fLambda_Daughter_Proton_ITS_nCluster(0),
  fLambda_Daughter_AntiPion_px(-999.0),
  fLambda_Daughter_AntiPion_py(-999.0),
  fLambda_Daughter_AntiPion_pz(-999.0),
  fLambda_Daughter_AntiPion_px_Generated(-999.0),
  fLambda_Daughter_AntiPion_py_Generated(-999.0),
  fLambda_Daughter_AntiPion_pz_Generated(-999.0),
  fLambda_Daughter_AntiPion_px_DecayVertex(-999.0),
  fLambda_Daughter_AntiPion_py_DecayVertex(-999.0),
  fLambda_Daughter_AntiPion_pz_DecayVertex(-999.0),
  fLambda_Daughter_AntiPion_pTPC(-999.0),
  fLambda_Daughter_AntiPion_Eta(-999.0),
  fLambda_Daughter_AntiPion_Phi(-999.0),
  fLambda_Daughter_AntiPion_TPC_Chi2(-999.0),
  fLambda_Daughter_AntiPion_TPC_dEdx(-999.0),
  fLambda_Daughter_AntiPion_TPC_dEdx_Sigma(-999.0),
  fLambda_Daughter_AntiPion_TOF_Mass2(-999.0),
  fLambda_Daughter_AntiPion_TOF_Mass2_Sigma(-999.0),
  fLambda_Daughter_AntiPion_ITS_dEdx(-999.0),
  fLambda_Daughter_AntiPion_ITS_dEdx_Sigma(-999.0),
  fLambda_Daughter_AntiPion_DCAxy(-999.0),
  fLambda_Daughter_AntiPion_DCAz(-999.0),
  fLambda_Daughter_AntiPion_TPC_nCrossedRows(0),
  fLambda_Daughter_AntiPion_TPC_nFindableCluster(0),
  fLambda_Daughter_AntiPion_TPC_nCluster(0),
  fLambda_Daughter_AntiPion_ITS_nCluster(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(-999.0),
  fDeuteron_py(-999.0),
  fDeuteron_pz(-999.0),
  fDeuteron_px_Generated(-999.0),
  fDeuteron_py_Generated(-999.0),
  fDeuteron_pz_Generated(-999.0),
  fDeuteron_pTPC(-999.0),
  fDeuteron_Eta(-999.0),
  fDeuteron_Phi(-999.0),
  fDeuteron_TPC_Chi2(-999.0),
  fDeuteron_TPC_dEdx(-999.0),
  fDeuteron_TPC_dEdx_Sigma(-999.0),
  fDeuteron_TOF_Mass2(-999.0),
  fDeuteron_TOF_Mass2_Sigma(-999.0),
  fDeuteron_ITS_dEdx(-999.0),
  fDeuteron_ITS_dEdx_Sigma(-999.0),
  fDeuteron_DCAxy(-999.0),
  fDeuteron_DCAz(-999.0),
  fDeuteron_Event_Centrality(-999.0),
  fDeuteron_Event_PrimaryVertexZ(-999.0),
  fDeuteron_Event_BField(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nFindableCluster(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_PDG(0),
  fDeuteron_MotherPDG(0),
  fDeuteron_FilterBit(0),
  fDeuteron_Event_Multiplicity(0),
  fDeuteron_Event_Identifier(0),
  fDeuteron_Event_RunNumber(0),
  fDeuteron_ITS_Layer0(0),
  fDeuteron_ITS_Layer1(0),
  fDeuteron_ITS_Layer2(0),
  fDeuteron_ITS_Layer3(0),
  fDeuteron_ITS_Layer4(0),
  fDeuteron_ITS_Layer5(0),
  fDeuteron_Event_IsFirstParticle(0),
  fSaveTree_AntiLambda(0),
  fAntiLambda_px(-999.0),
  fAntiLambda_py(-999.0),
  fAntiLambda_pz(-999.0),
  fAntiLambda_px_Generated(-999.0),
  fAntiLambda_py_Generated(-999.0),
  fAntiLambda_pz_Generated(-999.0),
  fAntiLambda_Eta(-999.0),
  fAntiLambda_Phi(-999.0),
  fAntiLambda_TransverseRadius(-999.0),
  fAntiLambda_CosinePointingAngle(-999.0),
  fAntiLambda_DCAv0ToPrimaryVertex(-999.0),
  fAntiLambda_DCAv0Daughters(-999.0),
  fAntiLambda_Alpha(-999.0),
  fAntiLambda_qT(-999.0),
  fAntiLambda_DecayLength(-999.0),
  fAntiLambda_PDG_Daughter1(0),
  fAntiLambda_PDG_Daughter2(0),
  fAntiLambda_PDG_v01(0),
  fAntiLambda_PDG_v02(0),
  fAntiLambda_PDG_Mother1(0),
  fAntiLambda_PDG_Mother2(0),
  fAntiLambda_SameV0(0),
  fAntiLambda_Event_Centrality(-999.0),
  fAntiLambda_Event_PrimaryVertexZ(-999.0),
  fAntiLambda_Event_BField(0),
  fAntiLambda_Event_Multiplicity(0),
  fAntiLambda_Event_Identifier(0),
  fAntiLambda_Event_RunNumber(0),
  fAntiLambda_Event_IsFirstParticle(0),
  fAntiLambda_Daughter_AntiProton_px(-999.0),
  fAntiLambda_Daughter_AntiProton_py(-999.0),
  fAntiLambda_Daughter_AntiProton_pz(-999.0),
  fAntiLambda_Daughter_AntiProton_px_Generated(-999.0),
  fAntiLambda_Daughter_AntiProton_py_Generated(-999.0),
  fAntiLambda_Daughter_AntiProton_pz_Generated(-999.0),
  fAntiLambda_Daughter_AntiProton_px_DecayVertex(-999.0),
  fAntiLambda_Daughter_AntiProton_py_DecayVertex(-999.0),
  fAntiLambda_Daughter_AntiProton_pz_DecayVertex(-999.0),
  fAntiLambda_Daughter_AntiProton_pTPC(-999.0),
  fAntiLambda_Daughter_AntiProton_Eta(-999.0),
  fAntiLambda_Daughter_AntiProton_Phi(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_Chi2(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2(-999.0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2_Sigma(-999.0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx(-999.0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_AntiProton_DCAxy(-999.0),
  fAntiLambda_Daughter_AntiProton_DCAz(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_AntiProton_TPC_nCluster(0),
  fAntiLambda_Daughter_AntiProton_ITS_nCluster(0),
  fAntiLambda_Daughter_Pion_px(-999.0),
  fAntiLambda_Daughter_Pion_py(-999.0),
  fAntiLambda_Daughter_Pion_pz(-999.0),
  fAntiLambda_Daughter_Pion_px_Generated(-999.0),
  fAntiLambda_Daughter_Pion_py_Generated(-999.0),
  fAntiLambda_Daughter_Pion_pz_Generated(-999.0),
  fAntiLambda_Daughter_Pion_px_DecayVertex(-999.0),
  fAntiLambda_Daughter_Pion_py_DecayVertex(-999.0),
  fAntiLambda_Daughter_Pion_pz_DecayVertex(-999.0),
  fAntiLambda_Daughter_Pion_pTPC(-999.0),
  fAntiLambda_Daughter_Pion_Eta(-999.0),
  fAntiLambda_Daughter_Pion_Phi(-999.0),
  fAntiLambda_Daughter_Pion_TPC_Chi2(-999.0),
  fAntiLambda_Daughter_Pion_TPC_dEdx(-999.0),
  fAntiLambda_Daughter_Pion_TPC_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_Pion_TOF_Mass2(-999.0),
  fAntiLambda_Daughter_Pion_TOF_Mass2_Sigma(-999.0),
  fAntiLambda_Daughter_Pion_ITS_dEdx(-999.0),
  fAntiLambda_Daughter_Pion_ITS_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_Pion_DCAxy(-999.0),
  fAntiLambda_Daughter_Pion_DCAz(-999.0),
  fAntiLambda_Daughter_Pion_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_Pion_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_Pion_TPC_nCluster(0),
  fAntiLambda_Daughter_Pion_ITS_nCluster(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(-999.0),
  fAntiDeuteron_py(-999.0),
  fAntiDeuteron_pz(-999.0),
  fAntiDeuteron_px_Generated(-999.0),
  fAntiDeuteron_py_Generated(-999.0),
  fAntiDeuteron_pz_Generated(-999.0),
  fAntiDeuteron_pTPC(-999.0),
  fAntiDeuteron_Eta(-999.0),
  fAntiDeuteron_Phi(-999.0),
  fAntiDeuteron_TPC_Chi2(-999.0),
  fAntiDeuteron_TPC_dEdx(-999.0),
  fAntiDeuteron_TPC_dEdx_Sigma(-999.0),
  fAntiDeuteron_TOF_Mass2(-999.0),
  fAntiDeuteron_TOF_Mass2_Sigma(-999.0),
  fAntiDeuteron_ITS_dEdx(-999.0),
  fAntiDeuteron_ITS_dEdx_Sigma(-999.0),
  fAntiDeuteron_DCAxy(-999.0),
  fAntiDeuteron_DCAz(-999.0),
  fAntiDeuteron_Event_Centrality(-999.0),
  fAntiDeuteron_Event_PrimaryVertexZ(-999.0),
  fAntiDeuteron_Event_BField(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nFindableCluster(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_PDG(0),
  fAntiDeuteron_MotherPDG(0),
  fAntiDeuteron_FilterBit(0),
  fAntiDeuteron_Event_Multiplicity(0),
  fAntiDeuteron_Event_Identifier(0),
  fAntiDeuteron_Event_RunNumber(0),
  fAntiDeuteron_ITS_Layer0(0),
  fAntiDeuteron_ITS_Layer1(0),
  fAntiDeuteron_ITS_Layer2(0),
  fAntiDeuteron_ITS_Layer3(0),
  fAntiDeuteron_ITS_Layer4(0),
  fAntiDeuteron_ITS_Layer5(0),
  fAntiDeuteron_Event_IsFirstParticle(0)
{


}



AliAnalysisTask_Ld_CreateTrees_PairsOnly::AliAnalysisTask_Ld_CreateTrees_PairsOnly(const char *name,Int_t CollisionSystem, Bool_t UseOpenCuts, Bool_t isMC, Bool_t SavePairsOnly) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(CollisionSystem),
  fUseOpenCuts(UseOpenCuts),
  fIsMC(isMC),
  fSavePairsOnly(SavePairsOnly),
  fSaveTree_Lambda(0),
  fLambda_px(-999.0),
  fLambda_py(-999.0),
  fLambda_pz(-999.0),
  fLambda_px_Generated(-999.0),
  fLambda_py_Generated(-999.0),
  fLambda_pz_Generated(-999.0),
  fLambda_Eta(-999.0),
  fLambda_Phi(-999.0),
  fLambda_TransverseRadius(-999.0),
  fLambda_CosinePointingAngle(-999.0),
  fLambda_DCAv0ToPrimaryVertex(-999.0),
  fLambda_DCAv0Daughters(-999.0),
  fLambda_Alpha(-999.0),
  fLambda_qT(-999.0),
  fLambda_DecayLength(-999.0),
  fLambda_PDG_Daughter1(0),
  fLambda_PDG_Daughter2(0),
  fLambda_PDG_v01(0),
  fLambda_PDG_v02(0),
  fLambda_PDG_Mother1(0),
  fLambda_PDG_Mother2(0),
  fLambda_SameV0(0),
  fLambda_Event_Centrality(-999.0),
  fLambda_Event_PrimaryVertexZ(-999.0),
  fLambda_Event_BField(0),
  fLambda_Event_Multiplicity(0),
  fLambda_Event_Identifier(0),
  fLambda_Event_RunNumber(0),
  fLambda_Event_IsFirstParticle(0),
  fLambda_Daughter_Proton_px(-999.0),
  fLambda_Daughter_Proton_py(-999.0),
  fLambda_Daughter_Proton_pz(-999.0),
  fLambda_Daughter_Proton_px_Generated(-999.0),
  fLambda_Daughter_Proton_py_Generated(-999.0),
  fLambda_Daughter_Proton_pz_Generated(-999.0),
  fLambda_Daughter_Proton_px_DecayVertex(-999.0),
  fLambda_Daughter_Proton_py_DecayVertex(-999.0),
  fLambda_Daughter_Proton_pz_DecayVertex(-999.0),
  fLambda_Daughter_Proton_pTPC(-999.0),
  fLambda_Daughter_Proton_Eta(-999.0),
  fLambda_Daughter_Proton_Phi(-999.0),
  fLambda_Daughter_Proton_TPC_Chi2(-999.0),
  fLambda_Daughter_Proton_TPC_dEdx(-999.0),
  fLambda_Daughter_Proton_TPC_dEdx_Sigma(-999.0),
  fLambda_Daughter_Proton_TOF_Mass2(-999.0),
  fLambda_Daughter_Proton_TOF_Mass2_Sigma(-999.0),
  fLambda_Daughter_Proton_ITS_dEdx(-999.0),
  fLambda_Daughter_Proton_ITS_dEdx_Sigma(-999.0),
  fLambda_Daughter_Proton_DCAxy(-999.0),
  fLambda_Daughter_Proton_DCAz(-999.0),
  fLambda_Daughter_Proton_TPC_nCrossedRows(0),
  fLambda_Daughter_Proton_TPC_nFindableCluster(0),
  fLambda_Daughter_Proton_TPC_nCluster(0),
  fLambda_Daughter_Proton_ITS_nCluster(0),
  fLambda_Daughter_AntiPion_px(-999.0),
  fLambda_Daughter_AntiPion_py(-999.0),
  fLambda_Daughter_AntiPion_pz(-999.0),
  fLambda_Daughter_AntiPion_px_Generated(-999.0),
  fLambda_Daughter_AntiPion_py_Generated(-999.0),
  fLambda_Daughter_AntiPion_pz_Generated(-999.0),
  fLambda_Daughter_AntiPion_px_DecayVertex(-999.0),
  fLambda_Daughter_AntiPion_py_DecayVertex(-999.0),
  fLambda_Daughter_AntiPion_pz_DecayVertex(-999.0),
  fLambda_Daughter_AntiPion_pTPC(-999.0),
  fLambda_Daughter_AntiPion_Eta(-999.0),
  fLambda_Daughter_AntiPion_Phi(-999.0),
  fLambda_Daughter_AntiPion_TPC_Chi2(-999.0),
  fLambda_Daughter_AntiPion_TPC_dEdx(-999.0),
  fLambda_Daughter_AntiPion_TPC_dEdx_Sigma(-999.0),
  fLambda_Daughter_AntiPion_TOF_Mass2(-999.0),
  fLambda_Daughter_AntiPion_TOF_Mass2_Sigma(-999.0),
  fLambda_Daughter_AntiPion_ITS_dEdx(-999.0),
  fLambda_Daughter_AntiPion_ITS_dEdx_Sigma(-999.0),
  fLambda_Daughter_AntiPion_DCAxy(-999.0),
  fLambda_Daughter_AntiPion_DCAz(-999.0),
  fLambda_Daughter_AntiPion_TPC_nCrossedRows(0),
  fLambda_Daughter_AntiPion_TPC_nFindableCluster(0),
  fLambda_Daughter_AntiPion_TPC_nCluster(0),
  fLambda_Daughter_AntiPion_ITS_nCluster(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(-999.0),
  fDeuteron_py(-999.0),
  fDeuteron_pz(-999.0),
  fDeuteron_px_Generated(-999.0),
  fDeuteron_py_Generated(-999.0),
  fDeuteron_pz_Generated(-999.0),
  fDeuteron_pTPC(-999.0),
  fDeuteron_Eta(-999.0),
  fDeuteron_Phi(-999.0),
  fDeuteron_TPC_Chi2(-999.0),
  fDeuteron_TPC_dEdx(-999.0),
  fDeuteron_TPC_dEdx_Sigma(-999.0),
  fDeuteron_TOF_Mass2(-999.0),
  fDeuteron_TOF_Mass2_Sigma(-999.0),
  fDeuteron_ITS_dEdx(-999.0),
  fDeuteron_ITS_dEdx_Sigma(-999.0),
  fDeuteron_DCAxy(-999.0),
  fDeuteron_DCAz(-999.0),
  fDeuteron_Event_Centrality(-999.0),
  fDeuteron_Event_PrimaryVertexZ(-999.0),
  fDeuteron_Event_BField(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nFindableCluster(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_PDG(0),
  fDeuteron_MotherPDG(0),
  fDeuteron_FilterBit(0),
  fDeuteron_Event_Multiplicity(0),
  fDeuteron_Event_Identifier(0),
  fDeuteron_Event_RunNumber(0),
  fDeuteron_ITS_Layer0(0),
  fDeuteron_ITS_Layer1(0),
  fDeuteron_ITS_Layer2(0),
  fDeuteron_ITS_Layer3(0),
  fDeuteron_ITS_Layer4(0),
  fDeuteron_ITS_Layer5(0),
  fDeuteron_Event_IsFirstParticle(0),
  fSaveTree_AntiLambda(0),
  fAntiLambda_px(-999.0),
  fAntiLambda_py(-999.0),
  fAntiLambda_pz(-999.0),
  fAntiLambda_px_Generated(-999.0),
  fAntiLambda_py_Generated(-999.0),
  fAntiLambda_pz_Generated(-999.0),
  fAntiLambda_Eta(-999.0),
  fAntiLambda_Phi(-999.0),
  fAntiLambda_TransverseRadius(-999.0),
  fAntiLambda_CosinePointingAngle(-999.0),
  fAntiLambda_DCAv0ToPrimaryVertex(-999.0),
  fAntiLambda_DCAv0Daughters(-999.0),
  fAntiLambda_Alpha(-999.0),
  fAntiLambda_qT(-999.0),
  fAntiLambda_DecayLength(-999.0),
  fAntiLambda_PDG_Daughter1(0),
  fAntiLambda_PDG_Daughter2(0),
  fAntiLambda_PDG_v01(0),
  fAntiLambda_PDG_v02(0),
  fAntiLambda_PDG_Mother1(0),
  fAntiLambda_PDG_Mother2(0),
  fAntiLambda_SameV0(0),
  fAntiLambda_Event_Centrality(-999.0),
  fAntiLambda_Event_PrimaryVertexZ(-999.0),
  fAntiLambda_Event_BField(0),
  fAntiLambda_Event_Multiplicity(0),
  fAntiLambda_Event_Identifier(0),
  fAntiLambda_Event_RunNumber(0),
  fAntiLambda_Event_IsFirstParticle(0),
  fAntiLambda_Daughter_AntiProton_px(-999.0),
  fAntiLambda_Daughter_AntiProton_py(-999.0),
  fAntiLambda_Daughter_AntiProton_pz(-999.0),
  fAntiLambda_Daughter_AntiProton_px_Generated(-999.0),
  fAntiLambda_Daughter_AntiProton_py_Generated(-999.0),
  fAntiLambda_Daughter_AntiProton_pz_Generated(-999.0),
  fAntiLambda_Daughter_AntiProton_px_DecayVertex(-999.0),
  fAntiLambda_Daughter_AntiProton_py_DecayVertex(-999.0),
  fAntiLambda_Daughter_AntiProton_pz_DecayVertex(-999.0),
  fAntiLambda_Daughter_AntiProton_pTPC(-999.0),
  fAntiLambda_Daughter_AntiProton_Eta(-999.0),
  fAntiLambda_Daughter_AntiProton_Phi(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_Chi2(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2(-999.0),
  fAntiLambda_Daughter_AntiProton_TOF_Mass2_Sigma(-999.0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx(-999.0),
  fAntiLambda_Daughter_AntiProton_ITS_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_AntiProton_DCAxy(-999.0),
  fAntiLambda_Daughter_AntiProton_DCAz(-999.0),
  fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_AntiProton_TPC_nCluster(0),
  fAntiLambda_Daughter_AntiProton_ITS_nCluster(0),
  fAntiLambda_Daughter_Pion_px(-999.0),
  fAntiLambda_Daughter_Pion_py(-999.0),
  fAntiLambda_Daughter_Pion_pz(-999.0),
  fAntiLambda_Daughter_Pion_px_Generated(-999.0),
  fAntiLambda_Daughter_Pion_py_Generated(-999.0),
  fAntiLambda_Daughter_Pion_pz_Generated(-999.0),
  fAntiLambda_Daughter_Pion_px_DecayVertex(-999.0),
  fAntiLambda_Daughter_Pion_py_DecayVertex(-999.0),
  fAntiLambda_Daughter_Pion_pz_DecayVertex(-999.0),
  fAntiLambda_Daughter_Pion_pTPC(-999.0),
  fAntiLambda_Daughter_Pion_Eta(-999.0),
  fAntiLambda_Daughter_Pion_Phi(-999.0),
  fAntiLambda_Daughter_Pion_TPC_Chi2(-999.0),
  fAntiLambda_Daughter_Pion_TPC_dEdx(-999.0),
  fAntiLambda_Daughter_Pion_TPC_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_Pion_TOF_Mass2(-999.0),
  fAntiLambda_Daughter_Pion_TOF_Mass2_Sigma(-999.0),
  fAntiLambda_Daughter_Pion_ITS_dEdx(-999.0),
  fAntiLambda_Daughter_Pion_ITS_dEdx_Sigma(-999.0),
  fAntiLambda_Daughter_Pion_DCAxy(-999.0),
  fAntiLambda_Daughter_Pion_DCAz(-999.0),
  fAntiLambda_Daughter_Pion_TPC_nCrossedRows(0),
  fAntiLambda_Daughter_Pion_TPC_nFindableCluster(0),
  fAntiLambda_Daughter_Pion_TPC_nCluster(0),
  fAntiLambda_Daughter_Pion_ITS_nCluster(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(-999.0),
  fAntiDeuteron_py(-999.0),
  fAntiDeuteron_pz(-999.0),
  fAntiDeuteron_px_Generated(-999.0),
  fAntiDeuteron_py_Generated(-999.0),
  fAntiDeuteron_pz_Generated(-999.0),
  fAntiDeuteron_pTPC(-999.0),
  fAntiDeuteron_Eta(-999.0),
  fAntiDeuteron_Phi(-999.0),
  fAntiDeuteron_TPC_Chi2(-999.0),
  fAntiDeuteron_TPC_dEdx(-999.0),
  fAntiDeuteron_TPC_dEdx_Sigma(-999.0),
  fAntiDeuteron_TOF_Mass2(-999.0),
  fAntiDeuteron_TOF_Mass2_Sigma(-999.0),
  fAntiDeuteron_ITS_dEdx(-999.0),
  fAntiDeuteron_ITS_dEdx_Sigma(-999.0),
  fAntiDeuteron_DCAxy(-999.0),
  fAntiDeuteron_DCAz(-999.0),
  fAntiDeuteron_Event_Centrality(-999.0),
  fAntiDeuteron_Event_PrimaryVertexZ(-999.0),
  fAntiDeuteron_Event_BField(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nFindableCluster(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_PDG(0),
  fAntiDeuteron_MotherPDG(0),
  fAntiDeuteron_FilterBit(0),
  fAntiDeuteron_Event_Multiplicity(0),
  fAntiDeuteron_Event_Identifier(0),
  fAntiDeuteron_Event_RunNumber(0),
  fAntiDeuteron_ITS_Layer0(0),
  fAntiDeuteron_ITS_Layer1(0),
  fAntiDeuteron_ITS_Layer2(0),
  fAntiDeuteron_ITS_Layer3(0),
  fAntiDeuteron_ITS_Layer4(0),
  fAntiDeuteron_ITS_Layer5(0),
  fAntiDeuteron_Event_IsFirstParticle(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());

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

}







void AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserCreateOutputObjects()
{


  fSaveTree_Lambda = new TTree("fSaveTree_Lambda","fSaveTree_Lambda");
  fSaveTree_Lambda->Branch("Lambda_px",&fLambda_px,"Lambda_px/F");
  fSaveTree_Lambda->Branch("Lambda_py",&fLambda_py,"Lambda_py/F");
  fSaveTree_Lambda->Branch("Lambda_pz",&fLambda_pz,"Lambda_pz/F");
  if(fIsMC == true){
  fSaveTree_Lambda->Branch("Lambda_px_Generated",&fLambda_px_Generated,"Lambda_px_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_py_Generated",&fLambda_py_Generated,"Lambda_py_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_pz_Generated",&fLambda_pz_Generated,"Lambda_pz_Generated/F");
  }
  fSaveTree_Lambda->Branch("Lambda_Eta",&fLambda_Eta,"Lambda_Eta/F");
  fSaveTree_Lambda->Branch("Lambda_Phi",&fLambda_Phi,"Lambda_Phi/F");
  fSaveTree_Lambda->Branch("Lambda_TransverseRadius",&fLambda_TransverseRadius,"Lambda_TransverseRadius/F");
  fSaveTree_Lambda->Branch("Lambda_CosinePointingAngle",&fLambda_CosinePointingAngle,"Lambda_CosinePointingAngle/F");
  fSaveTree_Lambda->Branch("Lambda_DCAv0ToPrimaryVertex",&fLambda_DCAv0ToPrimaryVertex,"Lambda_DCAv0ToPrimaryVertex/F");
  fSaveTree_Lambda->Branch("Lambda_DCAv0Daughters",&fLambda_DCAv0Daughters,"Lambda_DCAv0Daughters/F");
  fSaveTree_Lambda->Branch("Lambda_Alpha",&fLambda_Alpha,"Lambda_Alpha/F");
  fSaveTree_Lambda->Branch("Lambda_qT",&fLambda_qT,"Lambda_qT/F");
  fSaveTree_Lambda->Branch("Lambda_DecayLength",&fLambda_DecayLength,"Lambda_DecayLength/F");
  if(fIsMC == true){
  fSaveTree_Lambda->Branch("Lambda_PDG_Daughter1",&fLambda_PDG_Daughter1,"Lambda_PDG_Daughter1/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_Daughter2",&fLambda_PDG_Daughter2,"Lambda_PDG_Daughter2/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_v01",&fLambda_PDG_v01,"Lambda_PDG_v01/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_v02",&fLambda_PDG_v02,"Lambda_PDG_v02/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_Mother1",&fLambda_PDG_Mother1,"Lambda_PDG_Mother1/I");
  fSaveTree_Lambda->Branch("Lambda_PDG_Mother2",&fLambda_PDG_Mother2,"Lambda_PDG_Mother2/I");
  fSaveTree_Lambda->Branch("Lambda_SameV0",&fLambda_SameV0,"Lambda_SameV0/O");
  }
  fSaveTree_Lambda->Branch("Lambda_Event_Centrality",&fLambda_Event_Centrality,"Lambda_Event_Centrality/F");
  fSaveTree_Lambda->Branch("Lambda_Event_PrimaryVertexZ",&fLambda_Event_PrimaryVertexZ,"Lambda_Event_PrimaryVertexZ/F");
  fSaveTree_Lambda->Branch("Lambda_Event_BField",&fLambda_Event_BField,"Lambda_Event_BField/O");
  fSaveTree_Lambda->Branch("Lambda_Event_Multiplicity",&fLambda_Event_Multiplicity,"Lambda_Event_Multiplicity/I");
  fSaveTree_Lambda->Branch("Lambda_Event_Identifier",&fLambda_Event_Identifier,"Lambda_Event_Identifier/l");
  fSaveTree_Lambda->Branch("Lambda_Event_RunNumber",&fLambda_Event_RunNumber,"Lambda_Event_RunNumber/I");
  fSaveTree_Lambda->Branch("Lambda_Event_IsFirstParticle",&fLambda_Event_IsFirstParticle,"Lambda_Event_IsFirstParticle/O");

  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_px",&fLambda_Daughter_Proton_px,"Lambda_Daughter_Proton_px/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_py",&fLambda_Daughter_Proton_py,"Lambda_Daughter_Proton_py/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pz",&fLambda_Daughter_Proton_pz,"Lambda_Daughter_Proton_pz/F");
  if(fIsMC == true){
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_px_Generated",&fLambda_Daughter_Proton_px_Generated,"Lambda_Daughter_Proton_px_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_py_Generated",&fLambda_Daughter_Proton_py_Generated,"Lambda_Daughter_Proton_py_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pz_Generated",&fLambda_Daughter_Proton_pz_Generated,"Lambda_Daughter_Proton_pz_Generated/F");
  }
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_px_DecayVertex",&fLambda_Daughter_Proton_px_DecayVertex,"Lambda_Daughter_Proton_px_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_py_DecayVertex",&fLambda_Daughter_Proton_py_DecayVertex,"Lambda_Daughter_Proton_py_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pz_DecayVertex",&fLambda_Daughter_Proton_pz_DecayVertex,"Lambda_Daughter_Proton_pz_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_pTPC",&fLambda_Daughter_Proton_pTPC,"Lambda_Daughter_Proton_pTPC/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_Eta",&fLambda_Daughter_Proton_Eta,"Lambda_Daughter_Proton_Eta/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_Phi",&fLambda_Daughter_Proton_Phi,"Lambda_Daughter_Proton_Phi/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_Chi2",&fLambda_Daughter_Proton_TPC_Chi2,"Lambda_Daughter_Proton_TPC_Chi2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_dEdx",&fLambda_Daughter_Proton_TPC_dEdx,"Lambda_Daughter_Proton_TPC_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_dEdx_Sigma",&fLambda_Daughter_Proton_TPC_dEdx_Sigma,"Lambda_Daughter_Proton_TPC_dEdx_Sigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TOF_Mass2",&fLambda_Daughter_Proton_TOF_Mass2,"Lambda_Daughter_Proton_TOF_Mass2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TOF_Mass2_Sigma",&fLambda_Daughter_Proton_TOF_Mass2_Sigma,"Lambda_Daughter_Proton_TOF_Mass2_Sigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_dEdx",&fLambda_Daughter_Proton_ITS_dEdx,"Lambda_Daughter_Proton_ITS_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_dEdx_Sigma",&fLambda_Daughter_Proton_ITS_dEdx_Sigma,"Lambda_Daughter_Proton_ITS_dEdx_Sigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_DCAxy",&fLambda_Daughter_Proton_DCAxy,"Lambda_Daughter_Proton_DCAxy/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_DCAz",&fLambda_Daughter_Proton_DCAz,"Lambda_Daughter_Proton_DCAz/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nCrossedRows",&fLambda_Daughter_Proton_TPC_nCrossedRows,"Lambda_Daughter_Proton_TPC_nCrossedRows/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nFindableCluster",&fLambda_Daughter_Proton_TPC_nFindableCluster,"Lambda_Daughter_Proton_TPC_nFindableCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_TPC_nCluster",&fLambda_Daughter_Proton_TPC_nCluster,"Lambda_Daughter_Proton_TPC_nCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_Proton_ITS_nCluster",&fLambda_Daughter_Proton_ITS_nCluster,"Lambda_Daughter_Proton_ITS_nCluster/s");

  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_px",&fLambda_Daughter_AntiPion_px,"Lambda_Daughter_AntiPion_px/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_py",&fLambda_Daughter_AntiPion_py,"Lambda_Daughter_AntiPion_py/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz",&fLambda_Daughter_AntiPion_pz,"Lambda_Daughter_AntiPion_pz/F");
  if(fIsMC == true){
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_px_Generated",&fLambda_Daughter_AntiPion_px_Generated,"Lambda_Daughter_AntiPion_px_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_py_Generated",&fLambda_Daughter_AntiPion_py_Generated,"Lambda_Daughter_AntiPion_py_Generated/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz_Generated",&fLambda_Daughter_AntiPion_pz_Generated,"Lambda_Daughter_AntiPion_pz_Generated/F");
  }
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_px_DecayVertex",&fLambda_Daughter_AntiPion_px_DecayVertex,"Lambda_Daughter_AntiPion_px_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_py_DecayVertex",&fLambda_Daughter_AntiPion_py_DecayVertex,"Lambda_Daughter_AntiPion_py_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pz_DecayVertex",&fLambda_Daughter_AntiPion_pz_DecayVertex,"Lambda_Daughter_AntiPion_pz_DecayVertex/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_pTPC",&fLambda_Daughter_AntiPion_pTPC,"Lambda_Daughter_AntiPion_pTPC/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_Eta",&fLambda_Daughter_AntiPion_Eta,"Lambda_Daughter_AntiPion_Eta/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_Phi",&fLambda_Daughter_AntiPion_Phi,"Lambda_Daughter_AntiPion_Phi/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_Chi2",&fLambda_Daughter_AntiPion_TPC_Chi2,"Lambda_Daughter_AntiPion_TPC_Chi2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_dEdx",&fLambda_Daughter_AntiPion_TPC_dEdx,"Lambda_Daughter_AntiPion_TPC_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_dEdx_Sigma",&fLambda_Daughter_AntiPion_TPC_dEdx_Sigma,"Lambda_Daughter_AntiPion_TPC_dEdx_Sigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TOF_Mass2",&fLambda_Daughter_AntiPion_TOF_Mass2,"Lambda_Daughter_AntiPion_TOF_Mass2/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TOF_Mass2_Sigma",&fLambda_Daughter_AntiPion_TOF_Mass2_Sigma,"Lambda_Daughter_AntiPion_TOF_Mass2_Sigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_dEdx",&fLambda_Daughter_AntiPion_ITS_dEdx,"Lambda_Daughter_AntiPion_ITS_dEdx/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_dEdx_Sigma",&fLambda_Daughter_AntiPion_ITS_dEdx_Sigma,"Lambda_Daughter_AntiPion_ITS_dEdx_Sigma/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_DCAxy",&fLambda_Daughter_AntiPion_DCAxy,"Lambda_Daughter_AntiPion_DCAxy/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_DCAz",&fLambda_Daughter_AntiPion_DCAz,"Lambda_Daughter_AntiPion_DCAz/F");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nCrossedRows",&fLambda_Daughter_AntiPion_TPC_nCrossedRows,"Lambda_Daughter_AntiPion_TPC_nCrossedRows/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nFindableCluster",&fLambda_Daughter_AntiPion_TPC_nFindableCluster,"Lambda_Daughter_AntiPion_TPC_nFindableCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_TPC_nCluster",&fLambda_Daughter_AntiPion_TPC_nCluster,"Lambda_Daughter_AntiPion_TPC_nCluster/s");
  fSaveTree_Lambda->Branch("Lambda_Daughter_AntiPion_ITS_nCluster",&fLambda_Daughter_AntiPion_ITS_nCluster,"Lambda_Daughter_AntiPion_ITS_nCluster/s");









  fSaveTree_Deuteron = new TTree("fSaveTree_Deuteron","fSaveTree_Deuteron");
  fSaveTree_Deuteron->Branch("Deuteron_px",&fDeuteron_px,"Deuteron_px/F");
  fSaveTree_Deuteron->Branch("Deuteron_py",&fDeuteron_py,"Deuteron_py/F");
  fSaveTree_Deuteron->Branch("Deuteron_pz",&fDeuteron_pz,"Deuteron_pz/F");
  if(fIsMC == true){
  fSaveTree_Deuteron->Branch("Deuteron_px_Generated",&fDeuteron_px_Generated,"Deuteron_px_Generated/F");
  fSaveTree_Deuteron->Branch("Deuteron_py_Generated",&fDeuteron_py_Generated,"Deuteron_py_Generated/F");
  fSaveTree_Deuteron->Branch("Deuteron_pz_Generated",&fDeuteron_pz_Generated,"Deuteron_pz_Generated/F");
  }
  fSaveTree_Deuteron->Branch("Deuteron_pTPC",&fDeuteron_pTPC,"Deuteron_pTPC/F");
  fSaveTree_Deuteron->Branch("Deuteron_Eta",&fDeuteron_Eta,"Deuteron_Eta/F");
  fSaveTree_Deuteron->Branch("Deuteron_Phi",&fDeuteron_Phi,"Deuteron_Phi/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_Chi2",&fDeuteron_TPC_Chi2,"Deuteron_TPC_Chi2/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx",&fDeuteron_TPC_dEdx,"Deuteron_TPC_dEdx/F");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx_Sigma",&fDeuteron_TPC_dEdx_Sigma,"Deuteron_TPC_dEdx_Sigma/F");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2",&fDeuteron_TOF_Mass2,"Deuteron_TOF_Mass2/F");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2_Sigma",&fDeuteron_TOF_Mass2_Sigma,"Deuteron_TOF_Mass2_Sigma/F");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx",&fDeuteron_ITS_dEdx,"Deuteron_ITS_dEdx/F");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx_Sigma",&fDeuteron_ITS_dEdx_Sigma,"Deuteron_ITS_dEdx_Sigma/F");
  fSaveTree_Deuteron->Branch("Deuteron_DCAxy",&fDeuteron_DCAxy,"Deuteron_DCAxy/F");
  fSaveTree_Deuteron->Branch("Deuteron_DCAz",&fDeuteron_DCAz,"Deuteron_DCAz/F");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Centrality",&fDeuteron_Event_Centrality,"Deuteron_Event_Centrality/F");
  fSaveTree_Deuteron->Branch("Deuteron_Event_PrimaryVertexZ",&fDeuteron_Event_PrimaryVertexZ,"Deuteron_Event_PrimaryVertexZ/F");
  fSaveTree_Deuteron->Branch("Deuteron_Event_BField",&fDeuteron_Event_BField,"Deuteron_Event_BField/O");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&fDeuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nFindableCluster",&fDeuteron_TPC_nFindableCluster,"Deuteron_TPC_nFindableCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCluster",&fDeuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_nCluster",&fDeuteron_ITS_nCluster,"Deuteron_ITS_nCluster/s");
  if(fIsMC == true){
  fSaveTree_Deuteron->Branch("Deuteron_PDG",&fDeuteron_PDG,"Deuteron_PDG/I");
  fSaveTree_Deuteron->Branch("Deuteron_MotherPDG",&fDeuteron_MotherPDG,"Deuteron_MotherPDG/I");
  }
  fSaveTree_Deuteron->Branch("Deuteron_FilterBit",&fDeuteron_FilterBit,"Deuteron_FilterBit/O");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Multiplicity",&fDeuteron_Event_Multiplicity,"Deuteron_Event_Multiplicity/I");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Identifier",&fDeuteron_Event_Identifier,"Deuteron_Event_Identifier/l");
  fSaveTree_Deuteron->Branch("Deuteron_Event_RunNumber",&fDeuteron_Event_RunNumber,"Deuteron_Event_RunNumber/I");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_Layer0",&fDeuteron_ITS_Layer0,"Deuteron_ITS_Layer0/O");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_Layer1",&fDeuteron_ITS_Layer1,"Deuteron_ITS_Layer1/O");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_Layer2",&fDeuteron_ITS_Layer2,"Deuteron_ITS_Layer2/O");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_Layer3",&fDeuteron_ITS_Layer3,"Deuteron_ITS_Layer3/O");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_Layer4",&fDeuteron_ITS_Layer4,"Deuteron_ITS_Layer4/O");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_Layer5",&fDeuteron_ITS_Layer5,"Deuteron_ITS_Layer5/O");
  fSaveTree_Deuteron->Branch("Deuteron_Event_IsFirstParticle",&fDeuteron_Event_IsFirstParticle,"Deuteron_Event_IsFirstParticle/O");


  fSaveTree_AntiLambda = new TTree("fSaveTree_AntiLambda","fSaveTree_AntiLambda");
  fSaveTree_AntiLambda->Branch("AntiLambda_px",&fAntiLambda_px,"AntiLambda_px/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_py",&fAntiLambda_py,"AntiLambda_py/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_pz",&fAntiLambda_pz,"AntiLambda_pz/F");
  if(fIsMC == true){
  fSaveTree_AntiLambda->Branch("AntiLambda_px_Generated",&fAntiLambda_px_Generated,"AntiLambda_px_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_py_Generated",&fAntiLambda_py_Generated,"AntiLambda_py_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_pz_Generated",&fAntiLambda_pz_Generated,"AntiLambda_pz_Generated/F");
  }
  fSaveTree_AntiLambda->Branch("AntiLambda_Eta",&fAntiLambda_Eta,"AntiLambda_Eta/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Phi",&fAntiLambda_Phi,"AntiLambda_Phi/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_TransverseRadius",&fAntiLambda_TransverseRadius,"AntiLambda_TransverseRadius/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_CosinePointingAngle",&fAntiLambda_CosinePointingAngle,"AntiLambda_CosinePointingAngle/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_DCAv0ToPrimaryVertex",&fAntiLambda_DCAv0ToPrimaryVertex,"AntiLambda_DCAv0ToPrimaryVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_DCAv0Daughters",&fAntiLambda_DCAv0Daughters,"AntiLambda_DCAv0Daughters/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Alpha",&fAntiLambda_Alpha,"AntiLambda_Alpha/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_qT",&fAntiLambda_qT,"AntiLambda_qT/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_DecayLength",&fAntiLambda_DecayLength,"AntiLambda_DecayLength/F");
  if(fIsMC == true){
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Daughter1",&fAntiLambda_PDG_Daughter1,"AntiLambda_PDG_Daughter1/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Daughter2",&fAntiLambda_PDG_Daughter2,"AntiLambda_PDG_Daughter2/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_v01",&fAntiLambda_PDG_v01,"AntiLambda_PDG_v01/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_v02",&fAntiLambda_PDG_v02,"AntiLambda_PDG_v02/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Mother1",&fAntiLambda_PDG_Mother1,"AntiLambda_PDG_Mother1/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_PDG_Mother2",&fAntiLambda_PDG_Mother2,"AntiLambda_PDG_Mother2/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_SameV0",&fAntiLambda_SameV0,"AntiLambda_SameV0/O");
  }
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_Centrality",&fAntiLambda_Event_Centrality,"AntiLambda_Event_Centrality/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_PrimaryVertexZ",&fAntiLambda_Event_PrimaryVertexZ,"AntiLambda_Event_PrimaryVertexZ/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_BField",&fAntiLambda_Event_BField,"AntiLambda_Event_BField/O");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_Multiplicity",&fAntiLambda_Event_Multiplicity,"AntiLambda_Event_Multiplicity/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_Identifier",&fAntiLambda_Event_Identifier,"AntiLambda_Event_Identifier/l");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_RunNumber",&fAntiLambda_Event_RunNumber,"AntiLambda_Event_RunNumber/I");
  fSaveTree_AntiLambda->Branch("AntiLambda_Event_IsFirstParticle",&fAntiLambda_Event_IsFirstParticle,"AntiLambda_Event_IsFirstParticle/O");


  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px",&fAntiLambda_Daughter_AntiProton_px,"AntiLambda_Daughter_AntiProton_px/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py",&fAntiLambda_Daughter_AntiProton_py,"AntiLambda_Daughter_AntiProton_py/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz",&fAntiLambda_Daughter_AntiProton_pz,"AntiLambda_Daughter_AntiProton_pz/F");
  if(fIsMC == true){
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px_Generated",&fAntiLambda_Daughter_AntiProton_px_Generated,"AntiLambda_Daughter_AntiProton_px_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py_Generated",&fAntiLambda_Daughter_AntiProton_py_Generated,"AntiLambda_Daughter_AntiProton_py_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz_Generated",&fAntiLambda_Daughter_AntiProton_pz_Generated,"AntiLambda_Daughter_AntiProton_pz_Generated/F");
  }
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_px_DecayVertex",&fAntiLambda_Daughter_AntiProton_px_DecayVertex,"AntiLambda_Daughter_AntiProton_px_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_py_DecayVertex",&fAntiLambda_Daughter_AntiProton_py_DecayVertex,"AntiLambda_Daughter_AntiProton_py_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pz_DecayVertex",&fAntiLambda_Daughter_AntiProton_pz_DecayVertex,"AntiLambda_Daughter_AntiProton_pz_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_pTPC",&fAntiLambda_Daughter_AntiProton_pTPC,"AntiLambda_Daughter_AntiProton_pTPC/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_Eta",&fAntiLambda_Daughter_AntiProton_Eta,"AntiLambda_Daughter_AntiProton_Eta/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_Phi",&fAntiLambda_Daughter_AntiProton_Phi,"AntiLambda_Daughter_AntiProton_Phi/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_Chi2",&fAntiLambda_Daughter_AntiProton_TPC_Chi2,"AntiLambda_Daughter_AntiProton_TPC_Chi2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_dEdx",&fAntiLambda_Daughter_AntiProton_TPC_dEdx,"AntiLambda_Daughter_AntiProton_TPC_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_dEdx_Sigma",&fAntiLambda_Daughter_AntiProton_TPC_dEdx_Sigma,"AntiLambda_Daughter_AntiProton_TPC_dEdx_Sigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TOF_Mass2",&fAntiLambda_Daughter_AntiProton_TOF_Mass2,"AntiLambda_Daughter_AntiProton_TOF_Mass2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TOF_Mass2_Sigma",&fAntiLambda_Daughter_AntiProton_TOF_Mass2_Sigma,"AntiLambda_Daughter_AntiProton_TOF_Mass2_Sigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_dEdx",&fAntiLambda_Daughter_AntiProton_ITS_dEdx,"AntiLambda_Daughter_AntiProton_ITS_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_dEdx_Sigma",&fAntiLambda_Daughter_AntiProton_ITS_dEdx_Sigma,"AntiLambda_Daughter_AntiProton_ITS_dEdx_Sigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_DCAxy",&fAntiLambda_Daughter_AntiProton_DCAxy,"AntiLambda_Daughter_AntiProton_DCAxy/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_DCAz",&fAntiLambda_Daughter_AntiProton_DCAz,"AntiLambda_Daughter_AntiProton_DCAz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nCrossedRows",&fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows,"AntiLambda_Daughter_AntiProton_TPC_nCrossedRows/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nFindableCluster",&fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster,"AntiLambda_Daughter_AntiProton_TPC_nFindableCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_TPC_nCluster",&fAntiLambda_Daughter_AntiProton_TPC_nCluster,"AntiLambda_Daughter_AntiProton_TPC_nCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_AntiProton_ITS_nCluster",&fAntiLambda_Daughter_AntiProton_ITS_nCluster,"AntiLambda_Daughter_AntiProton_ITS_nCluster/s");

  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px",&fAntiLambda_Daughter_Pion_px,"AntiLambda_Daughter_Pion_px/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py",&fAntiLambda_Daughter_Pion_py,"AntiLambda_Daughter_Pion_py/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz",&fAntiLambda_Daughter_Pion_pz,"AntiLambda_Daughter_Pion_pz/F");
  if(fIsMC == true){
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px_Generated",&fAntiLambda_Daughter_Pion_px_Generated,"AntiLambda_Daughter_Pion_px_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py_Generated",&fAntiLambda_Daughter_Pion_py_Generated,"AntiLambda_Daughter_Pion_py_Generated/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz_Generated",&fAntiLambda_Daughter_Pion_pz_Generated,"AntiLambda_Daughter_Pion_pz_Generated/F");
  }
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_px_DecayVertex",&fAntiLambda_Daughter_Pion_px_DecayVertex,"AntiLambda_Daughter_Pion_px_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_py_DecayVertex",&fAntiLambda_Daughter_Pion_py_DecayVertex,"AntiLambda_Daughter_Pion_py_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pz_DecayVertex",&fAntiLambda_Daughter_Pion_pz_DecayVertex,"AntiLambda_Daughter_Pion_pz_DecayVertex/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_pTPC",&fAntiLambda_Daughter_Pion_pTPC,"AntiLambda_Daughter_Pion_pTPC/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_Eta",&fAntiLambda_Daughter_Pion_Eta,"AntiLambda_Daughter_Pion_Eta/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_Phi",&fAntiLambda_Daughter_Pion_Phi,"AntiLambda_Daughter_Pion_Phi/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_Chi2",&fAntiLambda_Daughter_Pion_TPC_Chi2,"AntiLambda_Daughter_Pion_TPC_Chi2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_dEdx",&fAntiLambda_Daughter_Pion_TPC_dEdx,"AntiLambda_Daughter_Pion_TPC_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_dEdx_Sigma",&fAntiLambda_Daughter_Pion_TPC_dEdx_Sigma,"AntiLambda_Daughter_Pion_TPC_dEdx_Sigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TOF_Mass2",&fAntiLambda_Daughter_Pion_TOF_Mass2,"AntiLambda_Daughter_Pion_TOF_Mass2/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TOF_Mass2_Sigma",&fAntiLambda_Daughter_Pion_TOF_Mass2_Sigma,"AntiLambda_Daughter_Pion_TOF_Mass2_Sigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_dEdx",&fAntiLambda_Daughter_Pion_ITS_dEdx,"AntiLambda_Daughter_Pion_ITS_dEdx/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_dEdx_Sigma",&fAntiLambda_Daughter_Pion_ITS_dEdx_Sigma,"AntiLambda_Daughter_Pion_ITS_dEdx_Sigma/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_DCAxy",&fAntiLambda_Daughter_Pion_DCAxy,"AntiLambda_Daughter_Pion_DCAxy/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_DCAz",&fAntiLambda_Daughter_Pion_DCAz,"AntiLambda_Daughter_Pion_DCAz/F");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nCrossedRows",&fAntiLambda_Daughter_Pion_TPC_nCrossedRows,"AntiLambda_Daughter_Pion_TPC_nCrossedRows/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nFindableCluster",&fAntiLambda_Daughter_Pion_TPC_nFindableCluster,"AntiLambda_Daughter_Pion_TPC_nFindableCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_TPC_nCluster",&fAntiLambda_Daughter_Pion_TPC_nCluster,"AntiLambda_Daughter_Pion_TPC_nCluster/s");
  fSaveTree_AntiLambda->Branch("AntiLambda_Daughter_Pion_ITS_nCluster",&fAntiLambda_Daughter_Pion_ITS_nCluster,"AntiLambda_Daughter_Pion_ITS_nCluster/s");



  fSaveTree_AntiDeuteron = new TTree("fSaveTree_AntiDeuteron","fSaveTree_AntiDeuteron");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_px",&fAntiDeuteron_px,"AntiDeuteron_px/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_py",&fAntiDeuteron_py,"AntiDeuteron_py/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pz",&fAntiDeuteron_pz,"AntiDeuteron_pz/F");
  if(fIsMC == true){
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_px_Generated",&fAntiDeuteron_px_Generated,"AntiDeuteron_px_Generated/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_py_Generated",&fAntiDeuteron_py_Generated,"AntiDeuteron_py_Generated/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pz_Generated",&fAntiDeuteron_pz_Generated,"AntiDeuteron_pz_Generated/F");
  }
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&fAntiDeuteron_pTPC,"AntiDeuteron_pTPC/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&fAntiDeuteron_Eta,"AntiDeuteron_Eta/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&fAntiDeuteron_Phi,"AntiDeuteron_Phi/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2",&fAntiDeuteron_TPC_Chi2,"AntiDeuteron_TPC_Chi2/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&fAntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_Sigma",&fAntiDeuteron_TPC_dEdx_Sigma,"AntiDeuteron_TPC_dEdx_Sigma/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&fAntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_Sigma",&fAntiDeuteron_TOF_Mass2_Sigma,"AntiDeuteron_TOF_Mass2_Sigma/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&fAntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_Sigma",&fAntiDeuteron_ITS_dEdx_Sigma,"AntiDeuteron_ITS_dEdx_Sigma/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&fAntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&fAntiDeuteron_DCAz,"AntiDeuteron_DCAz/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Centrality",&fAntiDeuteron_Event_Centrality,"AntiDeuteron_Event_Centrality/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_PrimaryVertexZ",&fAntiDeuteron_Event_PrimaryVertexZ,"AntiDeuteron_Event_PrimaryVertexZ/F");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_BField",&fAntiDeuteron_Event_BField,"AntiDeuteron_Event_BField/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&fAntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nFindableCluster",&fAntiDeuteron_TPC_nFindableCluster,"AntiDeuteron_TPC_nFindableCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&fAntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&fAntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/s");
  if(fIsMC == true){
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_PDG",&fAntiDeuteron_PDG,"AntiDeuteron_PDG/I");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_MotherPDG",&fAntiDeuteron_MotherPDG,"AntiDeuteron_MotherPDG/I");
  }
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_FilterBit",&fAntiDeuteron_FilterBit,"AntiDeuteron_FilterBit/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Multiplicity",&fAntiDeuteron_Event_Multiplicity,"AntiDeuteron_Event_Multiplicity/I");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&fAntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/l");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_RunNumber",&fAntiDeuteron_Event_RunNumber,"AntiDeuteron_Event_RunNumber/I");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_Layer0",&fAntiDeuteron_ITS_Layer0,"AntiDeuteron_ITS_Layer0/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_Layer1",&fAntiDeuteron_ITS_Layer1,"AntiDeuteron_ITS_Layer1/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_Layer2",&fAntiDeuteron_ITS_Layer2,"AntiDeuteron_ITS_Layer2/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_Layer3",&fAntiDeuteron_ITS_Layer3,"AntiDeuteron_ITS_Layer3/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_Layer4",&fAntiDeuteron_ITS_Layer4,"AntiDeuteron_ITS_Layer4/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_Layer5",&fAntiDeuteron_ITS_Layer5,"AntiDeuteron_ITS_Layer5/O");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_IsFirstParticle",&fAntiDeuteron_Event_IsFirstParticle,"AntiDeuteron_Event_IsFirstParticle/O");




  PostData(1,fSaveTree_Lambda);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiLambda);
  PostData(4,fSaveTree_AntiDeuteron);



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
  Bool_t DebugEventSelection = false;

  // define event cuts
  Float_t PrimaryVertexMaxZ = 10.0; // cm
  Float_t Centrality_min    = 0.0; // %
  Float_t Centrality_max    = 100.0; // %

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
  Int_t nTracks = fAODEvent->GetNumberOfTracks();
  if((nTracks == 0) || (nTracks < 0)) return;

  // use only events containing v0s
  Int_t nV0s = fAODEvent->GetNumberOfV0s();
  if((nV0s == 0) || (nV0s < 0)) return;

  // get primary vertex
  AliAODVertex *PrimaryVertex = fAODEvent->GetPrimaryVertex();
  if(!PrimaryVertex)::Warning("AliAnalsisTask_Ld_CreateTrees_PairsOnlyd::UserExec","No AliAODVertex object found!");
  Double_t PrimaryVertexPos[3] = {-999.0,-999.0,-999.0};
  PrimaryVertex->GetXYZ(PrimaryVertexPos);

  // apply cut on z-position of primary vertex
  Float_t PrimaryVertexZ = PrimaryVertexPos[2]; // cm
  if(TMath::IsNaN(PrimaryVertexZ)) return;
  if(TMath::Abs(PrimaryVertexZ) > PrimaryVertexMaxZ) return;

  // apply centrality cut
  Float_t Centrality = -999.0;
  AliMultSelection *MultSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
  Centrality = MultSelection->GetMultiplicityPercentile("V0M");
  if(TMath::IsNaN(Centrality)) return;

  if((fCollisionSystem == 1) || (fCollisionSystem == 2)){
    if((Centrality < Centrality_min) || (Centrality > Centrality_max)) return;
  }


  //combined reference multiplicity (tracklets + ITS + TPC) in |eta| < 0.8
  Int_t Multiplicity = fHeader->GetRefMultiplicityComb08();
  if((Multiplicity < 0) || (Multiplicity > 20000)) return;


  // get event information
  UInt_t PeriodNumber	    = fAODEvent->GetPeriodNumber();
  UInt_t OrbitNumber	    = fAODEvent->GetOrbitNumber();
  UShort_t BunchCrossNumber = fAODEvent->GetBunchCrossNumber();
  UInt_t TimeStamp	    = fAODEvent->GetTimeStamp();
  Int_t RunNumber	    = fAODEvent->GetRunNumber();
  Float_t BField	    = fAODEvent->GetMagneticField();

  if(RunNumber < 0) return;
  if(TMath::IsNaN(BField)) return;

  ULong64_t EventID = 0;
  ULong64_t Seed = 0;
  Int_t RandomNumber = 0;
  Int_t nTracksMC = 0;

  Bool_t IsPositiveBFieldPolarity = false;
  if(BField > 0.0) IsPositiveBFieldPolarity = true;

  // EventID -> https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventInfo
  if(fIsMC == false) EventID = ((ULong64_t)BunchCrossNumber) + ((ULong64_t)OrbitNumber*3564) + ((ULong64_t)PeriodNumber*16777215*3564);


  TRandom3 *RandomGenerator = new TRandom3();
  RandomGenerator->SetSeed(0);

  // EventID (MC)
  if(fIsMC == true){

    Seed = RandomGenerator->GetSeed();
    UInt_t Offset1 = 10000000;
    UInt_t Offset2 = 100000;
    RandomNumber = RandomGenerator->Integer(Offset2);

    Int_t IntegerPrimaryVertexZ = TMath::Abs((Int_t)(std::trunc(PrimaryVertexZ*1000000)));
    nTracksMC = fMCEvent->GetNumberOfTracks();
    EventID = ((((ULong64_t) RunNumber * Offset1) + (ULong64_t)IntegerPrimaryVertexZ) * Offset2 * 10) + (ULong64_t)RandomNumber;

  } // end of fIsMC == true








  // print event information
  if(DebugEventSelection)
  {

    cout << "" << endl;
    cout << "fCollisionSystem:\t\t" << fCollisionSystem << std::endl;
    cout << "PeriodNumber:\t\t\t" << PeriodNumber << endl;
    cout << "TimeStamp:\t\t\t" << TimeStamp << endl;
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
  // +++ Deuteron selection loop +++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Int_t nDeuteronsSelected = 0;
  Int_t nAntiDeuteronsSelected = 0;
  std::vector<Int_t> DeuteronVector;

  // Deuteron / Antideuteron loop
  for(Int_t track = 0; track < nTracks; track++){

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");

    Short_t Charge = 0;
    Charge = Track->Charge();
    if(!(TMath::Abs(Charge) == 1)) continue;

    Bool_t DeuteronIsMatter = false;
    if(Charge == +1) DeuteronIsMatter = true;
    if(Charge == -1) DeuteronIsMatter = false;

    Int_t ParticleSpecies = 0;
    if(DeuteronIsMatter == true) ParticleSpecies = 2;
    if(DeuteronIsMatter == false) ParticleSpecies = 4;

    Bool_t PassedDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,ParticleSpecies,RunNumber);
    if(PassedDeuteronCuts == false) continue;

    if(DeuteronIsMatter == true){

      DeuteronVector.push_back(track);
      nDeuteronsSelected++;

    } // end of isDeuteron

    if(DeuteronIsMatter == false){

      DeuteronVector.push_back(track);
      nAntiDeuteronsSelected++;

    } // end of is AntiDeuteron


  } // end of Deuteron / Antideuteron selection loop


  // if there are no Deuterons AND no Antideuterons get rid of this event
  if((fSavePairsOnly == true) && (nDeuteronsSelected == 0) && (nAntiDeuteronsSelected == 0)) return;










  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ Lambda selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TF1 *FitPionDCAxy_min = new TF1("FitPionDCAxy_min","pol2",0.0,6.0);
  FitPionDCAxy_min->FixParameter(0,0.0635603);
  FitPionDCAxy_min->FixParameter(1,-0.0593139);
  FitPionDCAxy_min->FixParameter(2,0.0171978);
  
  TF1 *FitPionDCAz_min = new TF1("FitPionDCAz_min","pol2",0.0,6.0);
  FitPionDCAz_min->FixParameter(0,0.071036);
  FitPionDCAz_min->FixParameter(1,-0.0441316);
  FitPionDCAz_min->FixParameter(2,0.00701049);

  Int_t nLambdasSelected = 0;
  Int_t nAntiLambdasSelected = 0;
  std::vector<Int_t> LambdaVector;
  std::vector<Int_t> AntilambdaVector;

  // Lambda / Antilambda loop
  for(Int_t V0 = 0; V0 < nV0s; V0++){ 

    AliAODv0 *v0 = dynamic_cast<AliAODv0*>(fAODEvent->GetV0(V0));
    if(!v0)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODv0 found");
    if(!v0) continue;

    AliAODTrack *PositiveTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    if(!PositiveTrack) continue;

    AliAODTrack *NegativeTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    if(!NegativeTrack) continue;

    Float_t PionDCAxy_min = -999.0;
    Float_t PionDCAz_min = -999.0;

    PionDCAxy_min = FitPionDCAxy_min->Eval(PositiveTrack->Pt());
    PionDCAz_min = FitPionDCAz_min->Eval(PositiveTrack->Pt());
    if(TMath::IsNaN(PionDCAxy_min)) continue;
    if(TMath::IsNaN(PionDCAz_min))  continue;
    Bool_t PassedProtonCuts = CheckProtonCuts(*PositiveTrack,*fPIDResponse,1,RunNumber);
    Bool_t PassedPionCuts = CheckPionCuts(*PositiveTrack,*fPIDResponse,5,RunNumber,PionDCAxy_min,PionDCAz_min);


    PionDCAxy_min = FitPionDCAxy_min->Eval(NegativeTrack->Pt());
    PionDCAz_min = FitPionDCAz_min->Eval(NegativeTrack->Pt());
    if(TMath::IsNaN(PionDCAxy_min)) continue;
    if(TMath::IsNaN(PionDCAz_min))  continue;
    Bool_t PassedAntiProtonCuts = CheckProtonCuts(*NegativeTrack,*fPIDResponse,3,RunNumber);
    Bool_t PassedAntiPionCuts = CheckPionCuts(*NegativeTrack,*fPIDResponse,6,RunNumber,PionDCAxy_min,PionDCAz_min); 

    if((PassedProtonCuts == false) && (PassedAntiProtonCuts == false)) continue;
    if((PassedPionCuts == false) && (PassedAntiPionCuts == false)) continue;
    if((PassedProtonCuts == true) && (PassedPionCuts == true)) continue;
    if((PassedAntiProtonCuts == true) && (PassedAntiPionCuts == true)) continue;

    Bool_t LambdaIsMatter = false;
    if((PassedProtonCuts == true) && (PassedAntiPionCuts == true)) LambdaIsMatter = true;
    if((PassedAntiProtonCuts == true) && (PassedPionCuts == true)) LambdaIsMatter = false;

    AliAODTrack *ProtonTrack = nullptr;
    AliAODTrack *PionTrack = nullptr;

    if(LambdaIsMatter == true){
    
      ProtonTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
      PionTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));

    }

    if(LambdaIsMatter == false){
    
      ProtonTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
      PionTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));

    }

    Short_t ProtonCharge = ProtonTrack->Charge();
    Short_t PionCharge = PionTrack->Charge();

    if(!(TMath::Abs(ProtonCharge) == 1)) continue;
    if(!(TMath::Abs(PionCharge) == 1)) continue;
    if(ProtonCharge == PionCharge) continue;


    AliAODVertex *DecayVertex = v0->GetSecondaryVtx();

    AliExternalTrackParam *ProtonTrackParam = new AliExternalTrackParam();
    ProtonTrackParam->CopyFromVTrack(ProtonTrack);
    Double_t Proton_DCA_new[2] = {-999.0,-999.0};
    Double_t Proton_CovarianceMatrix_new[3] = {-999.0,-999.0,-999.0};
    ProtonTrackParam->PropagateToDCA(DecayVertex,BField,10.0,Proton_DCA_new,Proton_CovarianceMatrix_new);
    Double_t ProtonMomentumPropagated[3] = {-999.0,-999.0,-999.0};
    ProtonTrackParam->GetPxPyPz(ProtonMomentumPropagated);

    AliExternalTrackParam *PionTrackParam = new AliExternalTrackParam();
    PionTrackParam->CopyFromVTrack(PionTrack);
    Double_t Pion_DCA_new[2] = {-999.0,-999.0};
    Double_t Pion_CovarianceMatrix_new[3] = {-999.0,-999.0};
    PionTrackParam->PropagateToDCA(DecayVertex,BField,10.0,Pion_DCA_new,Pion_CovarianceMatrix_new);
    Double_t PionMomentumPropagated[3] = {-999.0,-999.0,-999.0};
    PionTrackParam->GetPxPyPz(PionMomentumPropagated);


    Double_t MomentumDaughterPosX = 0.0;
    Double_t MomentumDaughterPosY = 0.0;
    Double_t MomentumDaughterPosZ = 0.0;
    Double_t MomentumDaughterNegX = 0.0;
    Double_t MomentumDaughterNegY = 0.0;
    Double_t MomentumDaughterNegZ = 0.0;
     
    if(LambdaIsMatter == true){

      MomentumDaughterPosX = ProtonMomentumPropagated[0];
      MomentumDaughterPosY = ProtonMomentumPropagated[1];
      MomentumDaughterPosZ = ProtonMomentumPropagated[2];
      MomentumDaughterNegX = PionMomentumPropagated[0];
      MomentumDaughterNegY = PionMomentumPropagated[1];
      MomentumDaughterNegZ = PionMomentumPropagated[2];

    } // end of LambdaIsMatter == true

     
    if(LambdaIsMatter == false){

      MomentumDaughterNegX = ProtonMomentumPropagated[0];
      MomentumDaughterNegY = ProtonMomentumPropagated[1];
      MomentumDaughterNegZ = ProtonMomentumPropagated[2];
      MomentumDaughterPosX = PionMomentumPropagated[0];
      MomentumDaughterPosY = PionMomentumPropagated[1];
      MomentumDaughterPosZ = PionMomentumPropagated[2];

    } // end of LambdaIsMatter == false


    Double_t MomentumV0X = v0->MomV0X();
    Double_t MomentumV0Y = v0->MomV0Y();
    Double_t MomentumV0Z = v0->MomV0Z();


    if(TMath::IsNaN(MomentumDaughterPosX)) continue;
    if(TMath::IsNaN(MomentumDaughterPosY)) continue;
    if(TMath::IsNaN(MomentumDaughterPosZ)) continue;
    if(TMath::IsNaN(MomentumDaughterNegX)) continue;
    if(TMath::IsNaN(MomentumDaughterNegY)) continue;
    if(TMath::IsNaN(MomentumDaughterNegZ)) continue;
    if(TMath::IsNaN(MomentumV0X)) continue;
    if(TMath::IsNaN(MomentumV0Y)) continue;
    if(TMath::IsNaN(MomentumV0Z)) continue;

    TVector3 *MomentumDaughterPositive = new TVector3();
    TVector3 *MomentumDaughterNegative = new TVector3();
    TVector3 *MomentumV0 = new TVector3();
        
    MomentumDaughterPositive->SetXYZ(MomentumDaughterPosX,MomentumDaughterPosY,MomentumDaughterPosZ);
    MomentumDaughterNegative->SetXYZ(MomentumDaughterNegX,MomentumDaughterNegY,MomentumDaughterNegZ);
    MomentumV0->SetXYZ(MomentumV0X,MomentumV0Y,MomentumV0Z);
  
    Double_t p_Longitudinal_Positive = MomentumDaughterPositive->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of positively charged daughter
    Double_t p_Longitudinal_Negative = MomentumDaughterNegative->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of negatively charged daughter
  
    Double_t Alpha = (p_Longitudinal_Positive - p_Longitudinal_Negative) / (p_Longitudinal_Positive + p_Longitudinal_Negative);
    Double_t qT = MomentumDaughterNegative->Perp(*MomentumV0);
    if(TMath::IsNaN(Alpha)) continue;
    if(TMath::IsNaN(qT)) continue;

    MomentumDaughterPositive->Delete();
    MomentumDaughterNegative->Delete();
    MomentumV0->Delete();             

    Float_t MassInvLambda = 0.0;
    Float_t MassInvWrongLambda = 0.0;
    Float_t MassInvKaonShort = 0.0;
    Float_t MassInvPhoton = 0.0;

    MassInvLambda = CalculateInvariantMassLambda(ProtonMomentumPropagated,PionMomentumPropagated,1);
    MassInvWrongLambda = CalculateInvariantMassLambda(ProtonMomentumPropagated,PionMomentumPropagated,2);
    MassInvKaonShort = CalculateInvariantMassLambda(ProtonMomentumPropagated,PionMomentumPropagated,3);
    MassInvPhoton = CalculateInvariantMassLambda(ProtonMomentumPropagated,PionMomentumPropagated,4);

    if(TMath::IsNaN(MassInvLambda)) continue;
    if(TMath::IsNaN(MassInvWrongLambda)) continue;
    if(TMath::IsNaN(MassInvKaonShort)) continue;
    if(TMath::IsNaN(MassInvPhoton)) continue;


    Bool_t PassedLambdaCuts = CheckLambdaCuts(*v0,PrimaryVertexPos,*fPIDResponse,LambdaIsMatter,RunNumber,MassInvLambda,MassInvWrongLambda,MassInvKaonShort,MassInvPhoton);
    if(PassedLambdaCuts == false) continue;
  

    if(LambdaIsMatter == true){

      LambdaVector.push_back(V0);
      nLambdasSelected++;

    } // end of isLambda

    if(LambdaIsMatter == false){

      AntilambdaVector.push_back(V0);
      nAntiLambdasSelected++;

    } // end of isAntiLambda


  } // end of Lambda / AntiLambda selection loop


  FitPionDCAxy_min->Delete();
  FitPionDCAz_min->Delete();


  // if there are no Lambdas AND no Antilambdas get rid of this event
  if((fSavePairsOnly == true) && (nLambdasSelected == 0) && (nAntiLambdasSelected == 0)) return;






  Bool_t SaveMatter = false;
  if((fSavePairsOnly == true) && (nDeuteronsSelected > 0) && (nLambdasSelected > 0)) SaveMatter = true;
  if((fSavePairsOnly == false) && ((nDeuteronsSelected > 0) || (nLambdasSelected > 0))) SaveMatter = true;

  Bool_t SaveAntiMatter = false;
  if((fSavePairsOnly == true) && (nAntiDeuteronsSelected > 0) && (nAntiLambdasSelected > 0)) SaveAntiMatter = true;
  if((fSavePairsOnly == false) && ((nAntiDeuteronsSelected > 0) || (nAntiLambdasSelected > 0))) SaveAntiMatter = true;


  if((SaveMatter == true) || (SaveAntiMatter == true)){


    Int_t nDeuteronsInVector = DeuteronVector.size();
    Bool_t IsFirstDeuteron = true;
    Bool_t IsFirstAntideuteron = true;

    // loop over Deuterons / Antideuterons
    for(Int_t iDeuteron = 0; iDeuteron < nDeuteronsInVector; iDeuteron++){

      Int_t track = DeuteronVector.at(iDeuteron);

      AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
      if(!Track)::Warning("AliAnalysisTask_Ld_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");


      AliMCParticle *MCParticle = 0x0;
      Int_t Label = TMath::Abs(Track->GetLabel());

      Float_t Generated_px = -999.0;
      Float_t Generated_py = -999.0;
      Float_t Generated_pz = -999.0;
      Int_t PDG = 0;
      Int_t MotherPDG = 0;

      if(fIsMC == true){

	MCParticle = (AliMCParticle*) fMCEvent->GetTrack(Label);
	PDG = MCParticle->PdgCode();
	if(MCParticle->IsPhysicalPrimary() == true) MotherPDG = 1;
	if(MCParticle->IsSecondaryFromMaterial() == true) MotherPDG = 2;
	if(MCParticle->IsSecondaryFromWeakDecay() == true){
      
	  Int_t LabelMother = TMath::Abs(MCParticle->GetMother());
	  AliMCParticle *MCParticleMother = (AliMCParticle*) fMCEvent->GetTrack(LabelMother);
	  MotherPDG = MCParticleMother->PdgCode();

	} // end of isSecondary

	Generated_px = MCParticle->Px();
	Generated_py = MCParticle->Py();
	Generated_pz = MCParticle->Pz();

      } // end of fIsMC == true


  
      Float_t xv[2] = {-999.0,-999.0};
      Float_t yv[3] = {-999.0,-999.0,-999.0};
      Track->GetImpactParameters(xv,yv);
      Float_t DCAxy = -999.0;
      Float_t DCAz = -999.0;
      DCAxy = xv[0];
      DCAz = xv[1];

      Int_t Charge = 0;
      Charge = Track->Charge();

      Bool_t DeuteronIsMatter = false;
      if(Charge > 0) DeuteronIsMatter = true;
      if(Charge < 0) DeuteronIsMatter = false;

      Int_t ParticleSpecies = 0;
      if(DeuteronIsMatter == true) ParticleSpecies = 2;
      if(DeuteronIsMatter == false) ParticleSpecies = 4;


      Float_t TOF_m2 = -999.0;
      Float_t TOF_m2_Sigma = -999.0;

      AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
      Bool_t TOFisOK = false;
      if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
      if(TOFisOK){

	TOF_m2 = CalculateMassSquareTOF(*Track);
	TOF_m2_Sigma = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,ParticleSpecies,RunNumber);

      } // end of TOFisOK




      Float_t ITS_dEdx = -999.0;
      Float_t ITS_dEdx_Sigma = -999.0;
      Int_t ITS_nCluster = 0;

      AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
      Bool_t ITSisOK = false;
      if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
      if(ITSisOK){

	ITS_dEdx = Track->GetITSsignal();
	ITS_dEdx_Sigma = CalculateSigmadEdxITS(*Track,ParticleSpecies,RunNumber);
	ITS_nCluster = Track->GetITSNcls();
	if((ITS_nCluster < 0) || (ITS_nCluster > 7)) ITS_nCluster = 0;

      } // end of ITSisOK


      if((DeuteronIsMatter == true) && (SaveMatter == true)){

	fDeuteron_px = Track->Px();
	fDeuteron_py = Track->Py();
	fDeuteron_pz = Track->Pz();
	fDeuteron_px_Generated = Generated_px;
	fDeuteron_py_Generated = Generated_py;
	fDeuteron_pz_Generated = Generated_pz;
	fDeuteron_pTPC = Track->GetTPCmomentum();
	fDeuteron_Eta = Track->Eta();
	fDeuteron_Phi = Track->Phi();
	fDeuteron_TPC_Chi2 = Track->GetTPCchi2();
	fDeuteron_TPC_dEdx = Track->GetTPCsignal();
	fDeuteron_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
	fDeuteron_TOF_Mass2 = (Float_t)TOF_m2;
	fDeuteron_TOF_Mass2_Sigma = (Float_t)TOF_m2_Sigma;
	fDeuteron_ITS_dEdx = (Float_t)ITS_dEdx;
	fDeuteron_ITS_dEdx_Sigma = (Float_t)ITS_dEdx_Sigma;
	fDeuteron_DCAxy = DCAxy;
	fDeuteron_DCAz = DCAz;
	fDeuteron_TPC_nCrossedRows = Track->GetTPCCrossedRows();
	fDeuteron_TPC_nFindableCluster = Track->GetTPCNclsF();
	fDeuteron_TPC_nCluster = Track->GetTPCNcls();
	fDeuteron_ITS_nCluster = (UShort_t)ITS_nCluster;
	fDeuteron_PDG = PDG;
	fDeuteron_MotherPDG = MotherPDG;
	fDeuteron_FilterBit = Track->TestFilterBit(BIT(0));
	fDeuteron_ITS_Layer0 = Track->HasPointOnITSLayer(0);
	fDeuteron_ITS_Layer1 = Track->HasPointOnITSLayer(1);
	fDeuteron_ITS_Layer2 = Track->HasPointOnITSLayer(2);
	fDeuteron_ITS_Layer3 = Track->HasPointOnITSLayer(3);
	fDeuteron_ITS_Layer4 = Track->HasPointOnITSLayer(4);
	fDeuteron_ITS_Layer5 = Track->HasPointOnITSLayer(5);
	fDeuteron_Event_Multiplicity = Multiplicity;
	fDeuteron_Event_Centrality = Centrality;
	fDeuteron_Event_PrimaryVertexZ = PrimaryVertexZ;
	fDeuteron_Event_BField = IsPositiveBFieldPolarity;
	fDeuteron_Event_Identifier = EventID;
	fDeuteron_Event_RunNumber = RunNumber;
	fDeuteron_Event_IsFirstParticle = IsFirstDeuteron;

	fSaveTree_Deuteron->Fill();

	IsFirstDeuteron = false;

      } // end of isDeuteron



      if((DeuteronIsMatter == false) && (SaveAntiMatter == true)){

	fAntiDeuteron_px = Track->Px();
	fAntiDeuteron_py = Track->Py();
	fAntiDeuteron_pz = Track->Pz();
	fAntiDeuteron_px_Generated = Generated_px;
	fAntiDeuteron_py_Generated = Generated_py;
	fAntiDeuteron_pz_Generated = Generated_pz;
	fAntiDeuteron_pTPC = Track->GetTPCmomentum();
	fAntiDeuteron_Eta = Track->Eta();
	fAntiDeuteron_Phi = Track->Phi();
	fAntiDeuteron_TPC_Chi2 = Track->GetTPCchi2();
	fAntiDeuteron_TPC_dEdx = Track->GetTPCsignal();
	fAntiDeuteron_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
	fAntiDeuteron_TOF_Mass2 = (Float_t)TOF_m2;
	fAntiDeuteron_TOF_Mass2_Sigma = (Float_t)TOF_m2_Sigma;
	fAntiDeuteron_ITS_dEdx = (Float_t)ITS_dEdx;
	fAntiDeuteron_ITS_dEdx_Sigma = (Float_t)ITS_dEdx_Sigma;
	fAntiDeuteron_DCAxy = DCAxy;
	fAntiDeuteron_DCAz = DCAz;
	fAntiDeuteron_TPC_nCrossedRows = Track->GetTPCCrossedRows();
	fAntiDeuteron_TPC_nFindableCluster = Track->GetTPCNclsF();
	fAntiDeuteron_TPC_nCluster = Track->GetTPCNcls();
	fAntiDeuteron_ITS_nCluster = (UShort_t)ITS_nCluster;
	fAntiDeuteron_PDG = PDG;
	fAntiDeuteron_MotherPDG = MotherPDG;
	fAntiDeuteron_FilterBit = Track->TestFilterBit(BIT(0));
	fAntiDeuteron_ITS_Layer0 = Track->HasPointOnITSLayer(0);
	fAntiDeuteron_ITS_Layer1 = Track->HasPointOnITSLayer(1);
	fAntiDeuteron_ITS_Layer2 = Track->HasPointOnITSLayer(2);
	fAntiDeuteron_ITS_Layer3 = Track->HasPointOnITSLayer(3);
	fAntiDeuteron_ITS_Layer4 = Track->HasPointOnITSLayer(4);
	fAntiDeuteron_ITS_Layer5 = Track->HasPointOnITSLayer(5);
	fAntiDeuteron_Event_Multiplicity = Multiplicity;
	fAntiDeuteron_Event_Centrality = Centrality;
	fAntiDeuteron_Event_PrimaryVertexZ = PrimaryVertexZ;
	fAntiDeuteron_Event_BField = IsPositiveBFieldPolarity;
	fAntiDeuteron_Event_Identifier = EventID;
	fAntiDeuteron_Event_RunNumber = RunNumber;
	fAntiDeuteron_Event_IsFirstParticle = IsFirstAntideuteron;

	fSaveTree_AntiDeuteron->Fill();

	IsFirstAntideuteron = false;

      } // end of isAntideuteron


    } // end of loop over Deuterons / Antideuterons







    Int_t nLambdasInVector = LambdaVector.size();
    Int_t nAntiLambdasInVector = AntilambdaVector.size();

    Bool_t IsFirstLambda = true;
    Bool_t IsFirstAntilambda = true;

    for(Int_t isMatter = 0; isMatter < 2; isMatter++){

      Int_t nParticlesInVector = 0;
      Bool_t LambdaIsMatter = false;

      if(isMatter == 0){
	
	if(nLambdasInVector == 0) continue;
	nParticlesInVector = nLambdasInVector;
	LambdaIsMatter = true;

      }

      if(isMatter == 1){

	if(nAntiLambdasInVector == 0) continue;
	nParticlesInVector = nAntiLambdasInVector;
	LambdaIsMatter = false;

      }


      // loop over Lambdas / Antilambdas
      for(Int_t iLambda = 0; iLambda < nParticlesInVector; iLambda++){

	Int_t V0 = 0;
	if(LambdaIsMatter == true)  V0 = LambdaVector.at(iLambda);
	if(LambdaIsMatter == false) V0 = AntilambdaVector.at(iLambda);

	AliAODv0 *v0 = dynamic_cast<AliAODv0*>(fAODEvent->GetV0(V0));
	AliAODTrack *ProtonTrack = nullptr;
	AliAODTrack *PionTrack = nullptr;

	Int_t ProtonSpecies = 0;
	Int_t PionSpecies = 0;

	if(LambdaIsMatter == true){

	  ProtonTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
	  PionTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
	  ProtonSpecies = 1;
	  PionSpecies = 6;

	}
	

	if(LambdaIsMatter == false){

	  ProtonTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
	  PionTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
	  ProtonSpecies = 3;
	  PionSpecies = 5;

	}

	


	AliAODVertex *DecayVertex = v0->GetSecondaryVtx();

	AliExternalTrackParam *ProtonTrackParam = new AliExternalTrackParam();
	ProtonTrackParam->CopyFromVTrack(ProtonTrack);
	Double_t Proton_DCA_new[2] = {-999.0,-999.0};
	Double_t Proton_CovarianceMatrix_new[3] = {-999.0,-999.0,-999.0};
	ProtonTrackParam->PropagateToDCA(DecayVertex,BField,10.0,Proton_DCA_new,Proton_CovarianceMatrix_new);
	Double_t ProtonMomentumPropagated[3] = {-999.0,-999.0,-999.0};
	ProtonTrackParam->GetPxPyPz(ProtonMomentumPropagated);

	AliExternalTrackParam *PionTrackParam = new AliExternalTrackParam();
	PionTrackParam->CopyFromVTrack(PionTrack);
	Double_t Pion_DCA_new[2] = {-999.0,-999.0};
	Double_t Pion_CovarianceMatrix_new[3] = {-999.0,-999.0};
	PionTrackParam->PropagateToDCA(DecayVertex,BField,10.0,Pion_DCA_new,Pion_CovarianceMatrix_new);
	Double_t PionMomentumPropagated[3] = {-999.0,-999.0,-999.0};
	PionTrackParam->GetPxPyPz(PionMomentumPropagated);


	Double_t MomentumDaughterPosX = 0.0;
	Double_t MomentumDaughterPosY = 0.0;
	Double_t MomentumDaughterPosZ = 0.0;
	Double_t MomentumDaughterNegX = 0.0;
	Double_t MomentumDaughterNegY = 0.0;
	Double_t MomentumDaughterNegZ = 0.0;
     
	if(LambdaIsMatter == true){

	  MomentumDaughterPosX = ProtonMomentumPropagated[0];
	  MomentumDaughterPosY = ProtonMomentumPropagated[1];
	  MomentumDaughterPosZ = ProtonMomentumPropagated[2];
	  MomentumDaughterNegX = PionMomentumPropagated[0];
	  MomentumDaughterNegY = PionMomentumPropagated[1];
	  MomentumDaughterNegZ = PionMomentumPropagated[2];

	} // end of LambdaIsMatter == true

     
	if(LambdaIsMatter == false){

	  MomentumDaughterNegX = ProtonMomentumPropagated[0];
	  MomentumDaughterNegY = ProtonMomentumPropagated[1];
	  MomentumDaughterNegZ = ProtonMomentumPropagated[2];
	  MomentumDaughterPosX = PionMomentumPropagated[0];
	  MomentumDaughterPosY = PionMomentumPropagated[1];
	  MomentumDaughterPosZ = PionMomentumPropagated[2];

	} // end of LambdaIsMatter == false

	Double_t MomentumV0X = v0->MomV0X();
	Double_t MomentumV0Y = v0->MomV0Y();
	Double_t MomentumV0Z = v0->MomV0Z();


	Float_t Generated_px = -999.0;
	Float_t Generated_py = -999.0;
	Float_t Generated_pz = -999.0;
	Float_t Generated_px_Daughter1 = -999.0;
	Float_t Generated_py_Daughter1 = -999.0;
	Float_t Generated_pz_Daughter1 = -999.0;
	Float_t Generated_px_Daughter2 = -999.0;
	Float_t Generated_py_Daughter2 = -999.0;
	Float_t Generated_pz_Daughter2 = -999.0;
	Int_t PDG_Daughter1 = 0;
	Int_t PDG_Daughter2 = 0;
	Int_t PDG_v01 = 0;
	Int_t PDG_v02 = 0;
	Int_t PDG_Mother1 = 0;
	Int_t PDG_Mother2 = 0;
	Bool_t SameV0 = false;

	Int_t Label_Daughter1 = TMath::Abs(ProtonTrack->GetLabel());
	Int_t Label_Daughter2 = TMath::Abs(PionTrack->GetLabel());


	if(fIsMC == true){

	  // get list of generated particles
	  TClonesArray *MCArray = dynamic_cast <TClonesArray*> (fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));

	  // Daughter1, v01, Mother1
	  AliMCParticle *MCParticleDaughter1 = (AliMCParticle*) MCArray->At(Label_Daughter1);
	  PDG_Daughter1 = MCParticleDaughter1->PdgCode();

	  Int_t Label_v01 = MCParticleDaughter1->GetMother();
	  AliMCParticle *MCv01 = (AliMCParticle*) MCArray->At(TMath::Abs(Label_v01));
	  PDG_v01 = MCv01->PdgCode();

	  Int_t Label_Mother1 = MCv01->GetMother();
	  AliMCParticle *MCParticleMother1 = (AliMCParticle*) MCArray->At(TMath::Abs(Label_Mother1));
	  PDG_Mother1 = MCParticleMother1->PdgCode();


	  // Daughter2, v02, Mother2
	  AliMCParticle *MCParticleDaughter2 = (AliMCParticle*) fMCEvent->GetTrack(Label_Daughter2);
	  PDG_Daughter2 = MCParticleDaughter2->PdgCode();

	  Int_t Label_v02 = MCParticleDaughter2->GetMother();
	  AliMCParticle *MCv02 = (AliMCParticle*) MCArray->At(TMath::Abs(Label_v02));
	  PDG_v02 = MCv02->PdgCode();

	  Int_t Label_Mother2 = MCv02->GetMother();
	  AliMCParticle *MCParticleMother2 = (AliMCParticle*) MCArray->At(TMath::Abs(Label_Mother2));
	  PDG_Mother2 = MCParticleMother2->PdgCode();



	  // check if daughters come from the same v0 
	  SameV0 = false;
	  if((Label_v01 == Label_v02) && (Label_v01 > 0) && (Label_v02 > 0)) SameV0 = true;


	  // where do the particles come from?
	  if(MCParticleDaughter1->IsPhysicalPrimary() == true) PDG_v01 = 1;
	  if(MCParticleDaughter1->IsSecondaryFromMaterial() == true) PDG_v01 = 2;

	  if(MCParticleDaughter2->IsPhysicalPrimary() == true) PDG_v02 = 1;
	  if(MCParticleDaughter2->IsSecondaryFromMaterial() == true) PDG_v02 = 2;

	  if(MCv01->IsPhysicalPrimary() == true) PDG_Mother1 = 1;
	  if(MCv01->IsSecondaryFromMaterial() == true) PDG_Mother1 = 2;

	  if(MCv02->IsPhysicalPrimary() == true) PDG_Mother2 = 1;
	  if(MCv02->IsSecondaryFromMaterial() == true) PDG_Mother2 = 2;

	  if(SameV0 == true){

	    // save generated momenta
	    Generated_px = MCv01->Px();
	    Generated_py = MCv01->Py();
	    Generated_pz = MCv01->Pz();

	  } // end of isSameV0

	  Generated_px_Daughter1 = MCParticleDaughter1->Px();
	  Generated_py_Daughter1 = MCParticleDaughter1->Py();
	  Generated_pz_Daughter1 = MCParticleDaughter1->Pz();

	  Generated_px_Daughter2 = MCParticleDaughter2->Px();
	  Generated_py_Daughter2 = MCParticleDaughter2->Py();
	  Generated_pz_Daughter2 = MCParticleDaughter2->Pz();


	} // end of fIsMC == true




  
	TVector3 *MomentumDaughterPositive = new TVector3();
	TVector3 *MomentumDaughterNegative = new TVector3();
	TVector3 *MomentumV0 = new TVector3();
	MomentumDaughterPositive->SetXYZ(MomentumDaughterPosX,MomentumDaughterPosY,MomentumDaughterPosZ);
	MomentumDaughterNegative->SetXYZ(MomentumDaughterNegX,MomentumDaughterNegY,MomentumDaughterNegZ);
	MomentumV0->SetXYZ(MomentumV0X,MomentumV0Y,MomentumV0Z);
  
	Double_t p_Longitudinal_Positive = MomentumDaughterPositive->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of positively charged daughter
	Double_t p_Longitudinal_Negative = MomentumDaughterNegative->Dot(*MomentumV0)/MomentumV0->Mag(); // longitudinal momentum of negatively charged daughter
  
	Double_t Alpha = (p_Longitudinal_Positive - p_Longitudinal_Negative) / (p_Longitudinal_Positive + p_Longitudinal_Negative);
	Double_t qT = MomentumDaughterNegative->Perp(*MomentumV0);


	MomentumDaughterPositive->Delete();
	MomentumDaughterNegative->Delete();
	MomentumV0->Delete(); 


	Float_t Proton_xv[2] = {-999.0,-999.0};
	Float_t Proton_yv[3] = {-999.0,-999.0,-999.0};
	ProtonTrack->GetImpactParameters(Proton_xv,Proton_yv);
	Float_t Proton_DCAxy = Proton_xv[0];
	Float_t Proton_DCAz = Proton_xv[1];



	Float_t Proton_TOF_m2 = -999.0;
	Float_t Proton_TOF_m2_Sigma = -999.0;

	AliPIDResponse::EDetPidStatus ProtonStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,ProtonTrack);
	Bool_t ProtonTOFisOK = false;
	if(ProtonStatusTOF == AliPIDResponse::kDetPidOk) ProtonTOFisOK = true;
    
	if(ProtonTOFisOK == true){

	  Proton_TOF_m2 = CalculateMassSquareTOF(*ProtonTrack);
	  Proton_TOF_m2_Sigma = CalculateSigmaMassSquareTOF(ProtonTrack->Pt(),Proton_TOF_m2,ProtonSpecies,RunNumber);

	}


	Float_t Proton_ITS_dEdx = -999.0;
	Float_t Proton_ITS_dEdx_Sigma = -999.0;
	Int_t Proton_ITS_nCluster = 0;

	AliPIDResponse::EDetPidStatus ProtonStatusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,ProtonTrack);
	Bool_t ProtonITSisOK = false;
	if(ProtonStatusITS == AliPIDResponse::kDetPidOk) ProtonITSisOK = true;
    
	if(ProtonITSisOK == true){

	  Proton_ITS_dEdx = ProtonTrack->GetITSsignal();
	  Proton_ITS_dEdx_Sigma = CalculateSigmadEdxITS(*ProtonTrack,ProtonSpecies,RunNumber);
	  Proton_ITS_nCluster = ProtonTrack->GetITSNcls();
	  if((Proton_ITS_nCluster < 0) || (Proton_ITS_nCluster > 7)) Proton_ITS_nCluster = 0;

	}



	Float_t Pion_xv[2] = {-999.0,-999.0};
	Float_t Pion_yv[3] = {-999.9,-999.0,-999.0};
	PionTrack->GetImpactParameters(Pion_xv,Pion_yv);
	Float_t Pion_DCAxy = Pion_xv[0];
	Float_t Pion_DCAz = Pion_xv[1];


	Float_t Pion_TOF_m2 = -999.0;
	Float_t Pion_TOF_m2_Sigma = -999.0;

	AliPIDResponse::EDetPidStatus PionStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,PionTrack);
	Bool_t PionTOFisOK = false;
	if(PionStatusTOF == AliPIDResponse::kDetPidOk) PionTOFisOK = true;
    
	if(PionTOFisOK == true){

	  Pion_TOF_m2 = CalculateMassSquareTOF(*PionTrack);
	  Pion_TOF_m2_Sigma = CalculateSigmaMassSquareTOF(PionTrack->Pt(),Pion_TOF_m2,PionSpecies,RunNumber);

	}


	Float_t Pion_ITS_dEdx = -999.0;
	Float_t Pion_ITS_dEdx_Sigma = -999.0;
	Int_t Pion_ITS_nCluster = 0;

	AliPIDResponse::EDetPidStatus PionStatusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,PionTrack);
	Bool_t PionITSisOK = false;
	if(PionStatusITS == AliPIDResponse::kDetPidOk) PionITSisOK = true;
    
	if(PionITSisOK == true){

	  Pion_ITS_dEdx = PionTrack->GetITSsignal();
	  Pion_ITS_dEdx_Sigma  = CalculateSigmadEdxITS(*PionTrack,PionSpecies,RunNumber);
	  Pion_ITS_nCluster     = PionTrack->GetITSNcls();
	  if((Pion_ITS_nCluster < 0) || (Pion_ITS_nCluster > 7)) Pion_ITS_nCluster = 0;

	}



	if((LambdaIsMatter == true) && (SaveMatter == true)){

	  fLambda_px = (Float_t)v0->MomV0X();
	  fLambda_py = (Float_t)v0->MomV0Y();
	  fLambda_pz = (Float_t)v0->MomV0Z();
	  fLambda_px_Generated = Generated_px;
	  fLambda_py_Generated = Generated_py;
	  fLambda_pz_Generated = Generated_pz;
	  fLambda_Eta = (Float_t)v0->PseudoRapV0();
	  fLambda_Phi = (Float_t)v0->Phi();
	  fLambda_TransverseRadius = (Float_t)v0->RadiusV0();
	  fLambda_CosinePointingAngle = (Float_t)TMath::Abs(v0->CosPointingAngle(PrimaryVertexPos));
	  fLambda_DCAv0ToPrimaryVertex = (Float_t)v0->DcaV0ToPrimVertex();
	  fLambda_DCAv0Daughters = (Float_t)v0->DcaV0Daughters();
	  fLambda_Alpha = (Float_t)Alpha;
	  fLambda_qT = (Float_t)qT;
	  fLambda_DecayLength = (Float_t)v0->DecayLengthV0(PrimaryVertexPos);
	  fLambda_PDG_Daughter1 = PDG_Daughter1;
	  fLambda_PDG_Daughter2 = PDG_Daughter2;
	  fLambda_PDG_v01 = PDG_v01;
	  fLambda_PDG_v02 = PDG_v02;
	  fLambda_PDG_Mother1 = PDG_Mother1;
	  fLambda_PDG_Mother2 = PDG_Mother2;
	  fLambda_SameV0 = SameV0;

	  fLambda_Daughter_Proton_px = ProtonTrack->Px();
	  fLambda_Daughter_Proton_py = ProtonTrack->Py();
	  fLambda_Daughter_Proton_pz = ProtonTrack->Pz();
	  fLambda_Daughter_Proton_px_Generated = Generated_px_Daughter1;
	  fLambda_Daughter_Proton_py_Generated = Generated_py_Daughter1;
	  fLambda_Daughter_Proton_pz_Generated = Generated_pz_Daughter1;
	  fLambda_Daughter_Proton_px_DecayVertex = MomentumDaughterPosX;
	  fLambda_Daughter_Proton_py_DecayVertex = MomentumDaughterPosY;
	  fLambda_Daughter_Proton_pz_DecayVertex = MomentumDaughterPosZ;
	  fLambda_Daughter_Proton_pTPC = ProtonTrack->GetTPCmomentum();
	  fLambda_Daughter_Proton_Eta = ProtonTrack->Eta();
	  fLambda_Daughter_Proton_Phi = ProtonTrack->Phi();
	  fLambda_Daughter_Proton_TPC_Chi2 = ProtonTrack->GetTPCchi2();
	  fLambda_Daughter_Proton_TPC_dEdx = ProtonTrack->GetTPCsignal();
	  fLambda_Daughter_Proton_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(ProtonTrack,AliPID::kProton);
	  fLambda_Daughter_Proton_TOF_Mass2 = (Float_t)Proton_TOF_m2;
	  fLambda_Daughter_Proton_TOF_Mass2_Sigma = (Float_t)Proton_TOF_m2_Sigma;
	  fLambda_Daughter_Proton_ITS_dEdx = (Float_t)Proton_ITS_dEdx;
	  fLambda_Daughter_Proton_ITS_dEdx_Sigma = (Float_t)Proton_ITS_dEdx_Sigma;
	  fLambda_Daughter_Proton_DCAxy = Proton_DCAxy;
	  fLambda_Daughter_Proton_DCAz = Proton_DCAz;
	  fLambda_Daughter_Proton_TPC_nCrossedRows = ProtonTrack->GetTPCCrossedRows();
	  fLambda_Daughter_Proton_TPC_nFindableCluster = ProtonTrack->GetTPCNclsF();
	  fLambda_Daughter_Proton_TPC_nCluster = ProtonTrack->GetTPCNcls();
	  fLambda_Daughter_Proton_ITS_nCluster = (UShort_t)Proton_ITS_nCluster;

	  fLambda_Daughter_AntiPion_px = PionTrack->Px();
	  fLambda_Daughter_AntiPion_py = PionTrack->Py();
	  fLambda_Daughter_AntiPion_pz = PionTrack->Pz();
	  fLambda_Daughter_AntiPion_px_Generated = Generated_px_Daughter2;
	  fLambda_Daughter_AntiPion_py_Generated = Generated_py_Daughter2;
	  fLambda_Daughter_AntiPion_pz_Generated = Generated_pz_Daughter2;
	  fLambda_Daughter_AntiPion_px_DecayVertex = MomentumDaughterNegX;
	  fLambda_Daughter_AntiPion_py_DecayVertex = MomentumDaughterNegY;
	  fLambda_Daughter_AntiPion_pz_DecayVertex = MomentumDaughterNegZ;
	  fLambda_Daughter_AntiPion_pTPC = PionTrack->GetTPCmomentum();
	  fLambda_Daughter_AntiPion_Eta = PionTrack->Eta();
	  fLambda_Daughter_AntiPion_Phi = PionTrack->Phi();
	  fLambda_Daughter_AntiPion_TPC_Chi2 = PionTrack->GetTPCchi2();
	  fLambda_Daughter_AntiPion_TPC_dEdx = PionTrack->GetTPCsignal();
	  fLambda_Daughter_AntiPion_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
	  fLambda_Daughter_AntiPion_TOF_Mass2 = (Float_t)Pion_TOF_m2;
	  fLambda_Daughter_AntiPion_TOF_Mass2_Sigma = (Float_t)Pion_TOF_m2_Sigma;
	  fLambda_Daughter_AntiPion_ITS_dEdx = (Float_t)Pion_ITS_dEdx;
	  fLambda_Daughter_AntiPion_ITS_dEdx_Sigma = (Float_t)Pion_ITS_dEdx_Sigma;
	  fLambda_Daughter_AntiPion_DCAxy = Pion_DCAxy;
	  fLambda_Daughter_AntiPion_DCAz = Pion_DCAz;
	  fLambda_Daughter_AntiPion_TPC_nCrossedRows = PionTrack->GetTPCCrossedRows();
	  fLambda_Daughter_AntiPion_TPC_nFindableCluster= PionTrack->GetTPCNclsF();
	  fLambda_Daughter_AntiPion_TPC_nCluster = PionTrack->GetTPCNcls();
	  fLambda_Daughter_AntiPion_ITS_nCluster = (UShort_t)Pion_ITS_nCluster;

	  fLambda_Event_Multiplicity = Multiplicity;
	  fLambda_Event_Centrality = Centrality;
	  fLambda_Event_PrimaryVertexZ = PrimaryVertexZ;
	  fLambda_Event_BField = IsPositiveBFieldPolarity;
	  fLambda_Event_Identifier = EventID;
	  fLambda_Event_RunNumber = RunNumber;
	  fLambda_Event_IsFirstParticle = IsFirstLambda;

	  fSaveTree_Lambda->Fill();

	  IsFirstLambda = false;

	} // end of SaveMatter



	if((LambdaIsMatter == false) && (SaveAntiMatter == true)){

	  fAntiLambda_px = (Float_t)v0->MomV0X();
	  fAntiLambda_py = (Float_t)v0->MomV0Y();
	  fAntiLambda_pz = (Float_t)v0->MomV0Z();
	  fAntiLambda_px_Generated = Generated_px;
	  fAntiLambda_py_Generated = Generated_py;
	  fAntiLambda_pz_Generated = Generated_pz;
	  fAntiLambda_Eta = (Float_t)v0->PseudoRapV0();
	  fAntiLambda_Phi = (Float_t)v0->Phi();
	  fAntiLambda_TransverseRadius = (Float_t)v0->RadiusV0();
	  fAntiLambda_CosinePointingAngle = (Float_t)TMath::Abs(v0->CosPointingAngle(PrimaryVertexPos));
	  fAntiLambda_DCAv0ToPrimaryVertex = (Float_t)v0->DcaV0ToPrimVertex();
	  fAntiLambda_DCAv0Daughters = (Float_t)v0->DcaV0Daughters();
	  fAntiLambda_Alpha = (Float_t)Alpha;
	  fAntiLambda_qT = (Float_t)qT;
	  fAntiLambda_DecayLength = (Float_t)v0->DecayLengthV0(PrimaryVertexPos);
	  fAntiLambda_PDG_Daughter1 = PDG_Daughter1;
	  fAntiLambda_PDG_Daughter2 = PDG_Daughter2;
	  fAntiLambda_PDG_v01 = PDG_v01;
	  fAntiLambda_PDG_v02 = PDG_v02;
	  fAntiLambda_PDG_Mother1 = PDG_Mother1;
	  fAntiLambda_PDG_Mother2 = PDG_Mother2;
	  fAntiLambda_SameV0 = SameV0;

	  fAntiLambda_Daughter_AntiProton_px = ProtonTrack->Px();
	  fAntiLambda_Daughter_AntiProton_py = ProtonTrack->Py();
	  fAntiLambda_Daughter_AntiProton_pz = ProtonTrack->Pz();
	  fAntiLambda_Daughter_AntiProton_px_Generated = Generated_px_Daughter1;
	  fAntiLambda_Daughter_AntiProton_py_Generated = Generated_py_Daughter1;
	  fAntiLambda_Daughter_AntiProton_pz_Generated = Generated_pz_Daughter1;
	  fAntiLambda_Daughter_AntiProton_px_DecayVertex = MomentumDaughterNegX;
	  fAntiLambda_Daughter_AntiProton_py_DecayVertex = MomentumDaughterNegY;
	  fAntiLambda_Daughter_AntiProton_pz_DecayVertex = MomentumDaughterNegZ;
	  fAntiLambda_Daughter_AntiProton_pTPC = ProtonTrack->GetTPCmomentum();
	  fAntiLambda_Daughter_AntiProton_Eta = ProtonTrack->Eta();
	  fAntiLambda_Daughter_AntiProton_Phi = ProtonTrack->Phi();
	  fAntiLambda_Daughter_AntiProton_TPC_Chi2 = ProtonTrack->GetTPCchi2();
	  fAntiLambda_Daughter_AntiProton_TPC_dEdx = ProtonTrack->GetTPCsignal();
	  fAntiLambda_Daughter_AntiProton_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(ProtonTrack,AliPID::kProton);
	  fAntiLambda_Daughter_AntiProton_TOF_Mass2 = (Float_t)Proton_TOF_m2;
	  fAntiLambda_Daughter_AntiProton_TOF_Mass2_Sigma = (Float_t)Proton_TOF_m2_Sigma;
	  fAntiLambda_Daughter_AntiProton_ITS_dEdx = (Float_t)Proton_ITS_dEdx;
	  fAntiLambda_Daughter_AntiProton_ITS_dEdx_Sigma = (Float_t)Proton_ITS_dEdx_Sigma;
	  fAntiLambda_Daughter_AntiProton_DCAxy = Proton_DCAxy;
	  fAntiLambda_Daughter_AntiProton_DCAz = Proton_DCAz;
	  fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows = ProtonTrack->GetTPCCrossedRows();
	  fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster = ProtonTrack->GetTPCNclsF();
	  fAntiLambda_Daughter_AntiProton_TPC_nCluster = ProtonTrack->GetTPCNcls();
	  fAntiLambda_Daughter_AntiProton_ITS_nCluster = (UShort_t)Proton_ITS_nCluster;

	  fAntiLambda_Daughter_Pion_px = PionTrack->Px();
	  fAntiLambda_Daughter_Pion_py = PionTrack->Py();
	  fAntiLambda_Daughter_Pion_pz = PionTrack->Pz();
	  fAntiLambda_Daughter_Pion_px_Generated = Generated_px_Daughter2;
	  fAntiLambda_Daughter_Pion_py_Generated = Generated_py_Daughter2;
	  fAntiLambda_Daughter_Pion_pz_Generated = Generated_pz_Daughter2;
	  fAntiLambda_Daughter_Pion_px_DecayVertex = MomentumDaughterPosX;
	  fAntiLambda_Daughter_Pion_py_DecayVertex = MomentumDaughterPosY;
	  fAntiLambda_Daughter_Pion_pz_DecayVertex = MomentumDaughterPosZ;
	  fAntiLambda_Daughter_Pion_pTPC = PionTrack->GetTPCmomentum();
	  fAntiLambda_Daughter_Pion_Eta = PionTrack->Eta();
	  fAntiLambda_Daughter_Pion_Phi = PionTrack->Phi();
	  fAntiLambda_Daughter_Pion_TPC_Chi2 = PionTrack->GetTPCchi2();
	  fAntiLambda_Daughter_Pion_TPC_dEdx = PionTrack->GetTPCsignal();
	  fAntiLambda_Daughter_Pion_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
	  fAntiLambda_Daughter_Pion_TOF_Mass2 = (Float_t)Pion_TOF_m2;
	  fAntiLambda_Daughter_Pion_TOF_Mass2_Sigma = (Float_t)Pion_TOF_m2_Sigma;
	  fAntiLambda_Daughter_Pion_ITS_dEdx = (Float_t)Pion_ITS_dEdx;
	  fAntiLambda_Daughter_Pion_ITS_dEdx_Sigma = (Float_t)Pion_ITS_dEdx_Sigma;
	  fAntiLambda_Daughter_Pion_DCAxy = Pion_DCAxy;
	  fAntiLambda_Daughter_Pion_DCAz = Pion_DCAz;
	  fAntiLambda_Daughter_Pion_TPC_nCrossedRows = PionTrack->GetTPCCrossedRows();
	  fAntiLambda_Daughter_Pion_TPC_nFindableCluster= PionTrack->GetTPCNclsF();
	  fAntiLambda_Daughter_Pion_TPC_nCluster = PionTrack->GetTPCNcls();
	  fAntiLambda_Daughter_Pion_ITS_nCluster = (UShort_t)Pion_ITS_nCluster;

	  fAntiLambda_Event_Multiplicity = Multiplicity;
	  fAntiLambda_Event_Centrality = Centrality;
	  fAntiLambda_Event_PrimaryVertexZ = PrimaryVertexZ;
	  fAntiLambda_Event_BField = IsPositiveBFieldPolarity;
	  fAntiLambda_Event_Identifier = EventID;
	  fAntiLambda_Event_RunNumber = RunNumber;
	  fAntiLambda_Event_IsFirstParticle = IsFirstAntilambda;

	  fSaveTree_AntiLambda->Fill();

	  IsFirstAntilambda = false;

	} // end of SaveMatter


      } // end of loop over Particles


    } // end of loop over matter / antimatter


  } // end of SaveMatter / SaveAntiMatter







  PostData(1,fSaveTree_Lambda);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiLambda);
  PostData(4,fSaveTree_AntiDeuteron);

} // end of UserExec















void AliAnalysisTask_Ld_CreateTrees_PairsOnly::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
Double_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateBetaTOF(AliAODTrack &track)
{

  Double_t length = track.GetIntegratedLength(); // cm
  if(TMath::IsNaN(length)) return -999.0;
  if(length <= 350.0) return -999.0;

  Double_t c = TMath::C(); // m/s
  Double_t end_time = track.GetTOFsignal(); // ps
  Double_t start_time = fPIDResponse->GetTOFResponse().GetStartTime(track.GetTPCmomentum()); // ps

  if(TMath::IsNaN(end_time)) return -999.0;
  if(TMath::IsNaN(start_time)) return -999.0;

  Double_t time = (end_time - start_time) * 1e-12; // ps -> s
  Double_t velocity = (length*0.01) / time; // m/s
  Double_t beta = velocity / c;

  return beta;

} // end of CalculateBetaTOF






// calculate the mass² TOF
Double_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateMassSquareTOF(AliAODTrack &track)
{

  Double_t mass2 = -999.0;

  Double_t p = track.P();
  Double_t beta = CalculateBetaTOF(track);

  if(TMath::IsNaN(p)) return mass2;
  if(TMath::IsNaN(beta)) return mass2;

  if(beta > 0.0){

    mass2 = (1/(beta*beta)-1) * (p*p);

  }

  return mass2;

} // end of CalculateMassSquareTOF








Double_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateSigmaMassSquareTOF(Double_t pT, Double_t massSq, Int_t ParticleSpecies, Int_t RunNumber)
{

  Double_t SigmaParticle = -999.0;
  if(TMath::IsNaN(massSq)) return SigmaParticle;
  if(massSq < -990.0) return SigmaParticle;


  Bool_t MetaLHC16 = false;
  Bool_t MetaLHC17 = false;
  Bool_t MetaLHC18 = false;
  Bool_t LHC18q = false;
  Bool_t LHC18r = false;

  if((RunNumber >= 252235) && (RunNumber <= 264347)) MetaLHC16 = true;
  if((RunNumber >= 270581) && (RunNumber <= 282704)) MetaLHC17 = true;
  if((RunNumber >= 285009) && (RunNumber <= 294925)) MetaLHC18 = true;
  if((RunNumber >= 295585) && (RunNumber <= 296623)) LHC18q = true;
  if((RunNumber >= 296690) && (RunNumber <= 297585)) LHC18r = true;


  Bool_t isProton	= false;
  Bool_t isDeuteron     = false;
  Bool_t isAntiProton   = false;
  Bool_t isAntiDeuteron = false;
  Bool_t isPion		= false;
  Bool_t isAntiPion     = false;

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




  Double_t mean = Mean->Eval(pT);
  Double_t sigma = Sigma->Eval(pT);

  Mean->Delete();
  Sigma->Delete();


  SigmaParticle = (massSq - mean)/(sigma);

  if(TMath::IsNaN(SigmaParticle)) return -999.0;

  return SigmaParticle;

} // end of CalculateSigmaMassSquareTOF










// apply track cuts for protons and antiprotons
bool AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, Int_t ParticleSpecies, Int_t RunNumber = 0)
{

  Bool_t PassedParticleCuts = false;



  // #######################################
  // ######## define particle cuts #########
  // #######################################

  Float_t Proton_pT_min = 0.0;
  Float_t Proton_pT_max = 999.0;
  Float_t Proton_eta_min = -0.9;
  Float_t Proton_eta_max = +0.9;
  Float_t Proton_p_max = 30.0;
  Float_t Proton_DCAxy_min = 0.01; // cm
  Float_t Proton_DCAxy_max = 999.0; // cm
  Float_t Proton_DCAz_min = 0.01; // cm
  Float_t Proton_DCAz_max = 999.0; // cm
  Float_t Proton_TPC_RatioRowsFindableCluster_min = 0.8;
  Float_t Proton_TPC_dEdx_Sigma_max = 4.0;
  Float_t Proton_TPC_Chi2perCluster_max = 5.0;
  Float_t Proton_TPC_Chi2perNDF_max = 4.0;
  Float_t Proton_TPC_Threshold = 1.0;
  Float_t Proton_TOF_Mass2_Sigma_max = 4.0;
  Float_t Proton_ITS_dEdx_Sigma_max = 4.0;
  Int_t Proton_TPC_nCluster_min = 50;
  Int_t Proton_TPC_nCrossedRows_min = 50;
  Int_t Proton_TPC_nSharedCluster_max = 0;
  Int_t Proton_ITS_nCluster_min = 0;
  Short_t Proton_Charge = 0;
  if(ParticleSpecies == 1) Proton_Charge = +1;
  if(ParticleSpecies == 3) Proton_Charge = -1;

  Bool_t UseTOF = true;
  Bool_t UseITS = true;


  // #######################################
  // ######## declare variables ############
  // #######################################

  Float_t p = -999.0;
  Float_t pT = -999.0;
  Float_t pTPC = -999.0;
  Float_t px = -999.0;
  Float_t py = -999.0;
  Float_t pz = -999.0;
  Float_t eta = -999.0;
  Float_t DCAxy = -999.0;
  Float_t DCAz = -999.0;
  Float_t TPC_dEdx = -999.0;
  Float_t TPC_dEdx_Sigma = -999.0;
  Float_t ITS_dEdx = -999.0;
  Float_t ITS_dEdx_Sigma = -999.0;
  Float_t TOF_Mass2 = -999.0;
  Float_t TOF_Mass2_Sigma = -999.0;
  Float_t TPC_RatioRowsFindableCluster = -999.0;
  Float_t TPC_Chi2 = -999.0;
  Float_t TPC_NDF = -999.0; 
  Float_t TPC_Chi2perCluster = -999.0;
  Float_t TPC_Chi2perNDF = -999.0;
  Float_t xv[2] = {-999.0,-999.0};
  Float_t yv[3] = {-999.0,-999.0,-999.0};

  Int_t TPC_nCluster = 0;
  Int_t TPC_nSharedCluster = 0;
  Int_t TPC_nCrossedRows = 0;
  Int_t TPC_nFindableCluster = 0;
  Int_t ITS_nCluster = 0;
  Short_t Charge = 0;

  Bool_t TPCisOK = false;
  Bool_t TOFisOK = false;
  Bool_t ITSisOK = false;




  // #######################################
  // ######## apply particle cuts ##########
  // #######################################

  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(statusTPC == AliPIDResponse::kDetPidOk) TPCisOK = true;
  if(TPCisOK == false) return PassedParticleCuts;

  p = Track.P();
  pT = Track.Pt();
  pTPC = Track.GetTPCmomentum();

  px = Track.Px();
  py = Track.Pz();
  pz = Track.Py();

  if(TMath::IsNaN(p)) return PassedParticleCuts;
  if(TMath::IsNaN(pT)) return PassedParticleCuts;
  if(TMath::IsNaN(pTPC)) return PassedParticleCuts;

  if(TMath::IsNaN(px)) return PassedParticleCuts;
  if(TMath::IsNaN(py)) return PassedParticleCuts;
  if(TMath::IsNaN(pz)) return PassedParticleCuts;

  if(TMath::Abs(px) > Proton_p_max) return PassedParticleCuts;
  if(TMath::Abs(py) > Proton_p_max) return PassedParticleCuts;
  if(TMath::Abs(pz) > Proton_p_max) return PassedParticleCuts;


  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Proton_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;

  // apply TPC Sigma cut
  TPC_dEdx_Sigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
  if(TMath::IsNaN(TPC_dEdx_Sigma)) return PassedParticleCuts;
  if(TMath::Abs(TPC_dEdx_Sigma) > Proton_TPC_dEdx_Sigma_max) return PassedParticleCuts;

  TPC_dEdx = Track.GetTPCsignal();
  if(TMath::IsNaN(TPC_dEdx)) return PassedParticleCuts;

  // apply charge cut
  Charge = Track.Charge();
  if(!(Proton_Charge == Charge)) return PassedParticleCuts;


  // get DCA information
  Track.GetImpactParameters(xv,yv);
  DCAxy = xv[0];
  DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;
  
  // apply DCAxy cut
  if((TMath::Abs(DCAxy) < Proton_DCAxy_min) || (TMath::Abs(DCAxy) > Proton_DCAxy_max)) return PassedParticleCuts;

  // apply DCAz cut
  if((TMath::Abs(DCAz) < Proton_DCAz_min) || (TMath::Abs(DCAz) > Proton_DCAz_max)) return PassedParticleCuts;

  // apply pT cut
  if(pT < Proton_pT_min || pT > Proton_pT_max) return PassedParticleCuts;


  // apply pseudo-rapidity cut
  eta = Track.Eta();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(eta < Proton_eta_min || eta > Proton_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  TPC_nCluster = Track.GetNcls(1);
  if((TMath::Abs(TPC_nCluster) > 170) || (TPC_nCluster < 0)) return PassedParticleCuts;
  if(TPC_nCluster < Proton_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  TPC_nCrossedRows = Track.GetTPCCrossedRows();
  if((TMath::Abs(TPC_nCrossedRows) > 170) || (TPC_nCluster < 0)) return PassedParticleCuts;
  if(TPC_nCrossedRows < Proton_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply shared cluster cut for TPC
  TPC_nSharedCluster = Track.GetTPCnclsS();
  if((TMath::Abs(TPC_nSharedCluster) > 10) || (TPC_nSharedCluster < 0)) return PassedParticleCuts;
  if(TPC_nSharedCluster > Proton_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  TPC_nFindableCluster = Track.GetTPCNclsF();
  if((TMath::Abs(TPC_nFindableCluster) > 170) || (TPC_nFindableCluster < 0)) return PassedParticleCuts;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((double)TPC_nCrossedRows / (double)TPC_nFindableCluster);
  if(TPC_RatioRowsFindableCluster < Proton_TPC_RatioRowsFindableCluster_min) return PassedParticleCuts;

  // receive TPC chi2
  TPC_Chi2 = Track.GetTPCchi2();
  if(TMath::IsNaN(TPC_Chi2)) return PassedParticleCuts;

  // compute TPC ndf
  if(TPC_nCluster > 5) TPC_NDF = (2*TPC_nCluster-5);

  // apply TPC chi2 per cluster cut
  TPC_Chi2perCluster = TPC_Chi2/((Float_t)TPC_nCluster);
  if(TPC_Chi2perCluster > Proton_TPC_Chi2perCluster_max) return PassedParticleCuts;

  // apply TPC chi2 per ndf cut
  TPC_Chi2perNDF = TPC_Chi2/TPC_NDF;
  if(TPC_Chi2perNDF > Proton_TPC_Chi2perNDF_max) return PassedParticleCuts; 




  // check if ITS information is available
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;


  if((ITSisOK == true) && (UseITS == true)){
    
    ITS_dEdx_Sigma = CalculateSigmadEdxITS(Track,ParticleSpecies,RunNumber);
    if(TMath::IsNaN(ITS_dEdx_Sigma)) return PassedParticleCuts;
    if(TMath::Abs(ITS_dEdx_Sigma) > Proton_ITS_dEdx_Sigma_max) return PassedParticleCuts;

    // apply ITS cluster cut
    ITS_nCluster = Track.GetITSNcls();
    if((TMath::Abs(ITS_nCluster) > 8) || (ITS_nCluster < 1)) return PassedParticleCuts;
    if(ITS_nCluster < Proton_ITS_nCluster_min) return PassedParticleCuts;

    ITS_dEdx = Track.GetITSsignal();
    if(TMath::IsNaN(ITS_dEdx)) return PassedParticleCuts;

  } // end of ITSisOK



  // apply TOF Mass2 cut
  if((TOFisOK == true) && (UseTOF == true)){

    TOF_Mass2 = CalculateMassSquareTOF(Track);
    if(TMath::IsNaN(TOF_Mass2)) return PassedParticleCuts;

    TOF_Mass2_Sigma = CalculateSigmaMassSquareTOF(pT,TOF_Mass2,ParticleSpecies,RunNumber);
    if(TMath::IsNaN(TOF_Mass2_Sigma)) return PassedParticleCuts;

    if(TMath::Abs(TOF_Mass2_Sigma) > Proton_TOF_Mass2_Sigma_max) return PassedParticleCuts;
 
  } // end of TOFisOK





  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckProtonCuts







// apply track cuts for Pions and antiPions
Bool_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckPionCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, Int_t ParticleSpecies = 0, Int_t RunNumber = 0, Float_t Pion_DCAxy_min = -999.0, Float_t Pion_DCAz_min = -999.0)
{

  Bool_t PassedParticleCuts = false;



  // #######################################
  // ######## define particle cuts #########
  // #######################################

  Float_t Pion_pT_min = 0.0;
  Float_t Pion_pT_max = 4.0;
  Float_t Pion_eta_min = -0.9;
  Float_t Pion_eta_max = +0.9;
  Float_t Pion_DCAxy_max = 999.0; // cm
  Float_t Pion_DCAz_max = 999.0; // cm
  Float_t Pion_p_max = 30.0;
  Float_t Pion_TPC_RatioRowsFindableCluster_min = 0.8;
  Float_t Pion_TPC_dEdx_Sigma_max = 4.0;
  Float_t Pion_TPC_Chi2perCluster_max = 5.0;
  Float_t Pion_TPC_Chi2perNDF_max = 4.0;
  Float_t Pion_TPC_Threshold = 4.0;
  Float_t Pion_TOF_Mass2_Sigma_max = 4.0;
  Float_t Pion_ITS_dEdx_Sigma_max = 4.0;
  Float_t Pion_ITS_nCluster_min = 0;
  Int_t Pion_TPC_nCluster_min = 50;
  Int_t Pion_TPC_nCrossedRows_min = 50;
  Int_t Pion_TPC_nSharedCluster_max = 0;
  Short_t Pion_Charge = 0;
  if(ParticleSpecies == 5) Pion_Charge = +1;
  if(ParticleSpecies == 6) Pion_Charge = -1;
  Bool_t UseTOF = true;
  Bool_t UseITS = true;



  // #######################################
  // ######## declare variables ############
  // #######################################

  Float_t p = -999.0;
  Float_t pT = -999.0;
  Float_t pTPC = -999.0;
  Float_t px = -999.0;
  Float_t py = -999.0;
  Float_t pz = -999.0;
  Float_t eta = -999.0;
  Float_t DCAxy = -999.0;
  Float_t DCAz = -999.0;
  Float_t TPC_dEdx = -999.0;
  Float_t TPC_dEdx_Sigma = -999.0;
  Float_t ITS_dEdx = -999.0;
  Float_t ITS_dEdx_Sigma = -999.0;
  Float_t TOF_Mass2 = -999.0;
  Float_t TOF_Mass2_Sigma = -999.0;
  Float_t TPC_RatioRowsFindableCluster = -999.0;
  Float_t TPC_Chi2 = -999.0;
  Float_t TPC_NDF = -999.0; 
  Float_t TPC_Chi2perCluster = -999.0;
  Float_t TPC_Chi2perNDF = -999.0;
  Float_t xv[2] = {-999.0,-999.0};
  Float_t yv[3] = {-999.0,-999.0,-999.0};

  Int_t TPC_nCluster = 0;
  Int_t TPC_nSharedCluster = 0;
  Int_t TPC_nCrossedRows = 0;
  Int_t TPC_nFindableCluster = 0;
  Int_t ITS_nCluster = 0;
  Short_t Charge = 0;

  Bool_t TPCisOK = false;
  Bool_t TOFisOK = false;
  Bool_t ITSisOK = false;




  // #######################################
  // ######## apply particle cuts ##########
  // #######################################

  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(statusTPC == AliPIDResponse::kDetPidOk) TPCisOK = true;
  if(TPCisOK == false) return PassedParticleCuts;

  p = Track.P();
  pT = Track.Pt();
  pTPC = Track.GetTPCmomentum();

  px = Track.Px();
  py = Track.Pz();
  pz = Track.Py();

  if(TMath::IsNaN(p)) return PassedParticleCuts;
  if(TMath::IsNaN(pT)) return PassedParticleCuts;
  if(TMath::IsNaN(pTPC)) return PassedParticleCuts;

  if(TMath::IsNaN(px)) return PassedParticleCuts;
  if(TMath::IsNaN(py)) return PassedParticleCuts;
  if(TMath::IsNaN(pz)) return PassedParticleCuts;

  if(TMath::Abs(px) > Pion_p_max) return PassedParticleCuts;
  if(TMath::Abs(py) > Pion_p_max) return PassedParticleCuts;
  if(TMath::Abs(pz) > Pion_p_max) return PassedParticleCuts;



  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Pion_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;

  // apply TPC Sigma cut
  TPC_dEdx_Sigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kPion);
  if(TMath::IsNaN(TPC_dEdx_Sigma)) return PassedParticleCuts;
  if(TMath::Abs(TPC_dEdx_Sigma) > Pion_TPC_dEdx_Sigma_max) return PassedParticleCuts;

  TPC_dEdx = Track.GetTPCsignal();
  if(TMath::IsNaN(TPC_dEdx)) return PassedParticleCuts;

  // get DCA information
  Track.GetImpactParameters(xv,yv);
  DCAxy = xv[0];
  DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;

  // apply charge cut 
  Charge = Track.Charge();
  if(!(Pion_Charge == Charge)) return PassedParticleCuts;
 
  // apply DCAxy cut
  if((TMath::Abs(DCAxy) < (2*Pion_DCAxy_min)) || (TMath::Abs(DCAxy) > Pion_DCAxy_max)) return PassedParticleCuts;

  // apply DCAz cut
  if((TMath::Abs(DCAz) < (2*Pion_DCAz_min)) || (TMath::Abs(DCAz) > Pion_DCAz_max)) return PassedParticleCuts;


  // apply pT cut
  if(pT < Pion_pT_min || pT > Pion_pT_max) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  eta = Track.Eta();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(eta < Pion_eta_min || eta > Pion_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  TPC_nCluster = Track.GetNcls(1);
  if((TMath::Abs(TPC_nCluster) > 170) || (TPC_nCluster < 0)) return PassedParticleCuts;
  if(TPC_nCluster < Pion_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  TPC_nCrossedRows = Track.GetTPCCrossedRows();
  if((TMath::Abs(TPC_nCrossedRows) > 170) || (TPC_nCluster < 0)) return PassedParticleCuts;
  if(TPC_nCrossedRows < Pion_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply shared cluster cut for TPC
  TPC_nSharedCluster = Track.GetTPCnclsS();
  if((TMath::Abs(TPC_nSharedCluster) > 10) || (TPC_nSharedCluster < 0)) return PassedParticleCuts;
  if(TPC_nSharedCluster > Pion_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  TPC_nFindableCluster = Track.GetTPCNclsF();
  if((TMath::Abs(TPC_nFindableCluster) > 170) || (TPC_nFindableCluster < 0)) return PassedParticleCuts;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((double)TPC_nCrossedRows / (double)TPC_nFindableCluster);
  if(TPC_RatioRowsFindableCluster < Pion_TPC_RatioRowsFindableCluster_min) return PassedParticleCuts;

  // receive TPC chi2
  TPC_Chi2 = Track.GetTPCchi2();
  if(TMath::IsNaN(TPC_Chi2)) return PassedParticleCuts;

  // compute TPC ndf
  if(TPC_nCluster > 5) TPC_NDF = (2*TPC_nCluster-5);

  // apply TPC chi2 per cluster cut
  TPC_Chi2perCluster = TPC_Chi2/((Float_t)TPC_nCluster);
  if(TPC_Chi2perCluster > Pion_TPC_Chi2perCluster_max) return PassedParticleCuts;

  // apply TPC chi2 per ndf cut
  TPC_Chi2perNDF = TPC_Chi2/TPC_NDF;
  if(TPC_Chi2perNDF > Pion_TPC_Chi2perNDF_max) return PassedParticleCuts; 




  // check if ITS information is available
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;

  if((ITSisOK == true) && (UseITS == true)){
  
    ITS_dEdx_Sigma = CalculateSigmadEdxITS(Track,ParticleSpecies,RunNumber);
    if(TMath::IsNaN(ITS_dEdx_Sigma)) return PassedParticleCuts;
    if(TMath::Abs(ITS_dEdx_Sigma) > Pion_ITS_dEdx_Sigma_max) return PassedParticleCuts;

    // apply ITS cluster cut
    ITS_nCluster = Track.GetITSNcls();
    if((TMath::Abs(ITS_nCluster) > 8) || (ITS_nCluster < 1)) return PassedParticleCuts;
    if(ITS_nCluster < Pion_ITS_nCluster_min) return PassedParticleCuts;

    ITS_dEdx = Track.GetITSsignal();
    if(TMath::IsNaN(ITS_dEdx)) return PassedParticleCuts;

  } // end of ITSisOK



  // apply TOF m2 cut   
  if((TOFisOK == true) && (UseTOF == true)){

    TOF_Mass2 = CalculateMassSquareTOF(Track);
    if(TMath::IsNaN(TOF_Mass2))  return PassedParticleCuts;

    TOF_Mass2_Sigma = CalculateSigmaMassSquareTOF(pT,TOF_Mass2,ParticleSpecies,RunNumber);
    if(TMath::IsNaN(TOF_Mass2_Sigma)) return PassedParticleCuts;

    if(TMath::Abs(TOF_Mass2_Sigma) > Pion_TOF_Mass2_Sigma_max) return PassedParticleCuts;
 
  } // end of TOFisOK





  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckPionCuts






































Bool_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckLambdaCuts(AliAODv0 &v0, Double_t PrimaryVertexPos[3], AliPIDResponse &fPIDResponse, Bool_t isMatter = true, Int_t RunNumber = 0, Float_t MassInvLambda = -999.0, Float_t MassInvWrongLambda = -999.0, Float_t MassInvKaonShort = -999.0, Float_t MassInvPhoton = -999.0)
{

  Bool_t PassedParticleCuts = false;



  // #######################################
  // ######## define particle cuts #########
  // #######################################

  Float_t Lambda_pT_min = 0.0; // GeV/c
  Float_t Lambda_pT_max = 3.0; // GeV/c
  Float_t Lambda_eta_min = -0.8;
  Float_t Lambda_eta_max = +0.8;
  Float_t Lambda_DecayRadius_min = 3.0; // cm
  Float_t Lambda_DecayRadius_max = 100.0; // cm
  Float_t Lambda_MassVariation_max = 0.01; // GeV/c²
  Float_t Lambda_DCAv0ToPrimaryVertex_min = 0.0; // cm
  Float_t Lambda_DCAv0ToPrimaryVertex_max = 1.0; // cm
  Float_t Lambda_p_max = 30.0;
  Float_t Lambda_CosinePointingAngle_min = 0.999;
  Float_t Lambda_DCAv0Daughters_min = 0.0;
  Float_t Lambda_DCAv0Daughters_max = 1.0;
  Float_t Lambda_MassPDG = 1.115683; // GeV/c²
  Bool_t Lambda_UseReconstructionOnTheFly = true;




  // #######################################
  // ######## declare variables ############
  // #######################################

  Bool_t IsReconstructedOnTheFly = false;
  Int_t nDaughters = 0;
  Int_t nProngs = 0;
  Short_t Charge = 0;
  Short_t Charge1 = 0;
  Short_t Charge2 = 0;
  Float_t eta = -999.0;
  Float_t DCAv0ToPrimaryVertex = -999.0;
  Float_t DCAv0Daughters = -999.0;
  Float_t DecayRadius = -999.0; 
  Float_t px = -999.0;
  Float_t py = -999.0;
  Float_t pz = -999.0;
  Float_t pT = -999.0;
  Float_t CosinePointingAngle = -999.0;



  // #######################################
  // ######## apply particle cuts ##########
  // #######################################

  IsReconstructedOnTheFly = v0.GetOnFlyStatus();
  if(!(IsReconstructedOnTheFly == Lambda_UseReconstructionOnTheFly)) return PassedParticleCuts;
 
  nDaughters = v0.GetNDaughters();
  if(!(nDaughters == 2)) return PassedParticleCuts;

  nProngs = v0.GetNProngs();
  if(!(nProngs == 2)) return PassedParticleCuts;

  Charge = v0.GetCharge();
  if(!(Charge == 0)) return PassedParticleCuts;

  Charge1 = v0.ChargeProng(0);
  Charge2 = v0.ChargeProng(1);
  if(!(TMath::Abs(Charge1) == 1)) return PassedParticleCuts;
  if(!(TMath::Abs(Charge2) == 1)) return PassedParticleCuts;
  if(Charge1 == Charge2) return PassedParticleCuts;
 
  eta = v0.PseudoRapV0();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if((eta < Lambda_eta_min) || (eta > Lambda_eta_max)) return PassedParticleCuts;


  DCAv0ToPrimaryVertex = v0.DcaV0ToPrimVertex();
  if(TMath::IsNaN(DCAv0ToPrimaryVertex)) return PassedParticleCuts;
  if((DCAv0ToPrimaryVertex < Lambda_DCAv0ToPrimaryVertex_min) || (DCAv0ToPrimaryVertex > Lambda_DCAv0ToPrimaryVertex_max)) return PassedParticleCuts;

  DecayRadius = v0.RadiusV0();
  if(TMath::IsNaN(DecayRadius)) return PassedParticleCuts;
  if((DecayRadius < Lambda_DecayRadius_min) || (DecayRadius > Lambda_DecayRadius_max)) return PassedParticleCuts;

  px = v0.MomV0X();
  py = v0.MomV0Y();
  pz = v0.MomV0Z();
  if(TMath::IsNaN(px)) return PassedParticleCuts;
  if(TMath::IsNaN(py)) return PassedParticleCuts;
  if(TMath::IsNaN(pz)) return PassedParticleCuts;

  if(TMath::Abs(px) > Lambda_p_max) return PassedParticleCuts;
  if(TMath::Abs(py) > Lambda_p_max) return PassedParticleCuts;
  if(TMath::Abs(pz) > Lambda_p_max) return PassedParticleCuts;


  pT = TMath::Sqrt((px * px) + (py * py));
  if((pT < Lambda_pT_min) || (pT > Lambda_pT_max)) return PassedParticleCuts; 

  CosinePointingAngle = v0.CosPointingAngle(PrimaryVertexPos);
  if(TMath::IsNaN(CosinePointingAngle)) return PassedParticleCuts;
  if(TMath::Abs(CosinePointingAngle) < Lambda_CosinePointingAngle_min) return PassedParticleCuts; 

  DCAv0Daughters = v0.DcaV0Daughters();
  if(TMath::IsNaN(DCAv0Daughters)) return PassedParticleCuts;
  if((DCAv0Daughters < Lambda_DCAv0Daughters_min) || (DCAv0Daughters > Lambda_DCAv0Daughters_max)) return PassedParticleCuts;



  if(fUseOpenCuts == true){

    if(MassInvLambda < 1.07 || MassInvLambda > 1.17) return PassedParticleCuts;
    if(MassInvWrongLambda < 1.4) return PassedParticleCuts;
    if(MassInvPhoton < 0.03) return PassedParticleCuts;
    if((MassInvKaonShort > 0.49) && (MassInvKaonShort < 0.51)) return PassedParticleCuts;

  }

  if(fUseOpenCuts == false){

    if(MassInvLambda < (Lambda_MassPDG-Lambda_MassVariation_max) || MassInvLambda > (Lambda_MassPDG+Lambda_MassVariation_max)) return PassedParticleCuts;
    if(MassInvWrongLambda < 1.4) return PassedParticleCuts;
    if(MassInvPhoton < 0.03) return PassedParticleCuts;
    if((MassInvKaonShort > 0.49) && (MassInvKaonShort < 0.51)) return PassedParticleCuts;

  }


  PassedParticleCuts = true;
  return PassedParticleCuts;

}







// apply track cuts for deuterons and antideuterons
Bool_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, Int_t ParticleSpecies = 0, Int_t RunNumber = 0)
{

  Bool_t PassedParticleCuts = false;


  // #######################################
  // ######## define particle cuts #########
  // #######################################

  Float_t Deuteron_pT_min = 0.0;
  Float_t Deuteron_pT_max = 3.0;
  Float_t Deuteron_eta_min = -0.8;
  Float_t Deuteron_eta_max = +0.8;
  Float_t Deuteron_DCAxy_max = 2.0; // cm
  Float_t Deuteron_DCAz_max = 0.1; // cm
  Float_t Deuteron_p_max = 30.0;
  Float_t Deuteron_TPC_RatioRowsFindableCluster_min = 0.8;
  Float_t Deuteron_TPC_dEdx_Sigma_max = 4.0;
  Float_t Deuteron_TPC_Chi2perCluster_max = 5.0;
  Float_t Deuteron_TPC_Chi2perNDF_max = 4.0;
  Int_t Deuteron_TPC_nCluster_min = 70;
  Int_t Deuteron_TPC_nCrossedRows_min = 60;
  Int_t Deuteron_TPC_nSharedCluster_max = 0;
  Float_t Deuteron_TPC_Threshold = 1.4;
  Float_t Deuteron_TOF_Mass2_Sigma_max = 4.0;
  Float_t Deuteron_ITS_dEdx_Sigma_max = 4.0;
  Int_t Deuteron_ITS_nCluster_min = 0;
  Bool_t UseTOF = true;
  Bool_t UseITS = true;

  Bool_t Extend_pT_range = true;
  Float_t pT_ExtensionThreshold = 1.6;
  Float_t Pion_TPC_dEdx_Sigma_max     = 3.0;
  Float_t Kaon_TPC_dEdx_Sigma_max     = 3.0;
  Float_t Proton_TPC_dEdx_Sigma_max   = 3.0;
  Float_t Electron_TPC_dEdx_Sigma_max = 3.0;
  Float_t Muon_TPC_dEdx_Sigma_max     = 3.0;



  // #######################################
  // ######## declare variables ############
  // #######################################

  Float_t p = -999.0;
  Float_t pT = -999.0;
  Float_t pTPC = -999.0;
  Float_t px = -999.0;
  Float_t py = -999.0;
  Float_t pz = -999.0;
  Float_t eta = -999.0;
  Float_t DCAxy = -999.0;
  Float_t DCAz = -999.0;
  Float_t TPC_dEdx = -999.0;
  Float_t TPC_dEdx_Sigma = -999.0;
  Float_t ITS_dEdx = -999.0;
  Float_t ITS_dEdx_Sigma = -999.0;
  Float_t TOF_Mass2 = -999.0;
  Float_t TOF_Mass2_Sigma = -999.0;
  Float_t TPC_RatioRowsFindableCluster = -999.0;
  Float_t TPC_Chi2 = -999.0;
  Float_t TPC_NDF = -999.0; 
  Float_t TPC_Chi2perCluster = -999.0;
  Float_t TPC_Chi2perNDF = -999.0;
  Float_t xv[2] = {-999.0,-999.0};
  Float_t yv[3] = {-999.0,-999.0,-999.0};

  Int_t TPC_nCluster = 0;
  Int_t TPC_nSharedCluster = 0;
  Int_t TPC_nCrossedRows = 0;
  Int_t TPC_nFindableCluster = 0;
  Int_t ITS_nCluster = 0;
  Short_t Charge = 0;

  Bool_t TPCisOK = false;
  Bool_t TOFisOK = false;
  Bool_t ITSisOK = false;

  Float_t TPC_dEdx_Sigma_Pion = -999.0;
  Float_t TPC_dEdx_Sigma_Kaon = -999.0;
  Float_t TPC_dEdx_Sigma_Proton = -999.0;
  Float_t TPC_dEdx_Sigma_Electron = -999.0;
  Float_t TPC_dEdx_Sigma_Muon = -999.0;



  // #######################################
  // ######## apply particle cuts ##########
  // #######################################

  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(statusTPC == AliPIDResponse::kDetPidOk) TPCisOK = true;
  if(TPCisOK == false) return PassedParticleCuts;

  p = Track.P();
  pT = Track.Pt();
  pTPC = Track.GetTPCmomentum();

  px = Track.Px();
  py = Track.Pz();
  pz = Track.Py();

  if(TMath::IsNaN(p)) return PassedParticleCuts;
  if(TMath::IsNaN(pT)) return PassedParticleCuts;
  if(TMath::IsNaN(pTPC)) return PassedParticleCuts;

  if(TMath::IsNaN(px)) return PassedParticleCuts;
  if(TMath::IsNaN(py)) return PassedParticleCuts;
  if(TMath::IsNaN(pz)) return PassedParticleCuts;

  if(TMath::Abs(px) > Deuteron_p_max) return PassedParticleCuts;
  if(TMath::Abs(py) > Deuteron_p_max) return PassedParticleCuts;
  if(TMath::Abs(pz) > Deuteron_p_max) return PassedParticleCuts;



  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Deuteron_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;


  // apply TPC Sigma cut
  TPC_dEdx_Sigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
  if(TMath::IsNaN(TPC_dEdx_Sigma)) return PassedParticleCuts;

  if(fUseOpenCuts == false){
    if(TMath::Abs(TPC_dEdx_Sigma) > Deuteron_TPC_dEdx_Sigma_max) return PassedParticleCuts;
  }

  if(fUseOpenCuts == true){
    if(TMath::Abs(TPC_dEdx_Sigma) > 7.0) return PassedParticleCuts;
  }

  TPC_dEdx = Track.GetTPCsignal();
  if(TMath::IsNaN(TPC_dEdx)) return PassedParticleCuts;

  // get DCA information
  Track.GetImpactParameters(xv,yv);
  DCAxy = xv[0];
  DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;
  
  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Deuteron_DCAxy_max) return PassedParticleCuts;

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Deuteron_DCAz_max) return PassedParticleCuts;

  // apply pT cut
  if(pT < Deuteron_pT_min || pT > Deuteron_pT_max) return PassedParticleCuts;

  Charge = Track.Charge();
  if(!(TMath::Abs(Charge) == 1)) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  eta = Track.Eta();
  if(TMath::IsNaN(eta)) return PassedParticleCuts;
  if(eta < Deuteron_eta_min || eta > Deuteron_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  TPC_nCluster = Track.GetNcls(1);
  if((TMath::Abs(TPC_nCluster) > 170) || (TPC_nCluster < 0)) return PassedParticleCuts;
  if(TPC_nCluster < Deuteron_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  TPC_nCrossedRows = Track.GetTPCCrossedRows();
  if((TMath::Abs(TPC_nCrossedRows) > 170) || (TPC_nCrossedRows < 0)) return PassedParticleCuts;
  if(TPC_nCrossedRows < Deuteron_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply zero shared cluster cut for TPC
  TPC_nSharedCluster = Track.GetTPCnclsS();
  if((TMath::Abs(TPC_nSharedCluster) > 10) || (TPC_nSharedCluster < 0)) return PassedParticleCuts;
  if(TPC_nSharedCluster > Deuteron_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  TPC_nFindableCluster = Track.GetTPCNclsF();
  if((TMath::Abs(TPC_nFindableCluster) > 170) || (TPC_nFindableCluster < 0)) return PassedParticleCuts;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((Float_t)TPC_nCrossedRows / (Float_t)TPC_nFindableCluster);
  if(TPC_RatioRowsFindableCluster < Deuteron_TPC_RatioRowsFindableCluster_min) return PassedParticleCuts;

  // receive TPC chi2
  TPC_Chi2 = Track.GetTPCchi2();
  if(TMath::IsNaN(TPC_Chi2)) return PassedParticleCuts;

  // compute TPC ndf
  if(TPC_nCluster > 5) TPC_NDF = (2*TPC_nCluster-5);

  // apply TPC chi2 per cluster cut
  TPC_Chi2perCluster = TPC_Chi2/((Float_t)TPC_nCluster);
  if(TPC_Chi2perCluster > Deuteron_TPC_Chi2perCluster_max) return PassedParticleCuts;

  // apply TPC chi2 per ndf cut
  TPC_Chi2perNDF = TPC_Chi2/TPC_NDF;
  if(TPC_Chi2perNDF > Deuteron_TPC_Chi2perNDF_max) return PassedParticleCuts; 


  if(Extend_pT_range == true){

    // cut out dEdx band of other particles above pT = 1.6 GeV/c²

    TPC_dEdx_Sigma_Pion = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kPion);
    if(TMath::IsNaN(TPC_dEdx_Sigma_Pion)) return PassedParticleCuts;
    if((pT >= pT_ExtensionThreshold) && (TMath::Abs(TPC_dEdx_Sigma_Pion) < Pion_TPC_dEdx_Sigma_max)) return PassedParticleCuts;

    TPC_dEdx_Sigma_Kaon = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kKaon);
    if(TMath::IsNaN(TPC_dEdx_Sigma_Kaon)) return PassedParticleCuts;
    if((pT >= pT_ExtensionThreshold) && (TMath::Abs(TPC_dEdx_Sigma_Kaon) < Kaon_TPC_dEdx_Sigma_max)) return PassedParticleCuts;

    TPC_dEdx_Sigma_Proton = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
    if(TMath::IsNaN(TPC_dEdx_Sigma_Proton)) return PassedParticleCuts;
    if((pT >= pT_ExtensionThreshold) && (TMath::Abs(TPC_dEdx_Sigma_Proton) < Proton_TPC_dEdx_Sigma_max)) return PassedParticleCuts;

    TPC_dEdx_Sigma_Electron = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kElectron);
    if(TMath::IsNaN(TPC_dEdx_Sigma_Electron)) return PassedParticleCuts;
    if((pT >= pT_ExtensionThreshold) && (TMath::Abs(TPC_dEdx_Sigma_Electron) < Electron_TPC_dEdx_Sigma_max)) return PassedParticleCuts;

    TPC_dEdx_Sigma_Muon = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kMuon);
    if(TMath::IsNaN(TPC_dEdx_Sigma_Muon)) return PassedParticleCuts;
    if((pT >= pT_ExtensionThreshold) && (TMath::Abs(TPC_dEdx_Sigma_Muon) < Muon_TPC_dEdx_Sigma_max)) return PassedParticleCuts;

  } // end of Extend_pT_range


  // check if ITS information is available
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;



  if((ITSisOK == true) && (UseITS == true)){
    
    ITS_dEdx_Sigma = CalculateSigmadEdxITS(Track,ParticleSpecies,RunNumber);
    if(TMath::IsNaN(ITS_dEdx_Sigma)) return PassedParticleCuts;
    if(TMath::Abs(ITS_dEdx_Sigma) > Deuteron_ITS_dEdx_Sigma_max) return PassedParticleCuts;

    // apply ITS cluster cut
    ITS_nCluster = Track.GetITSNcls();
    if((TMath::Abs(ITS_nCluster) > 8) || (ITS_nCluster < 1)) return PassedParticleCuts;
    if(ITS_nCluster < Deuteron_ITS_nCluster_min) return PassedParticleCuts;

    ITS_dEdx = Track.GetITSsignal();
    if(TMath::IsNaN(ITS_dEdx)) return PassedParticleCuts;

  } // end of ITSisOK






  if((TOFisOK == true) && (UseTOF == true)){

    TOF_Mass2 = CalculateMassSquareTOF(Track);
    if(TMath::IsNaN(TOF_Mass2)) return PassedParticleCuts;

    TOF_Mass2_Sigma = CalculateSigmaMassSquareTOF(pT,TOF_Mass2,ParticleSpecies,RunNumber);
    if(TMath::IsNaN(TOF_Mass2_Sigma)) return PassedParticleCuts;

    if(TMath::Abs(TOF_Mass2_Sigma) > Deuteron_TOF_Mass2_Sigma_max) return PassedParticleCuts;

  } // end of TOFisOK



  PassedParticleCuts = true;
  return PassedParticleCuts;


} // end of CheckDeuteronCuts
















Double_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateSigmadEdxITS(AliAODTrack &Track, Int_t ParticleSpecies, Int_t RunNumber){

  Double_t SigmaParticle = -999.0;
  Double_t SignalITS = Track.GetITSsignal();
  if(TMath::IsNaN(SignalITS)) return SigmaParticle;

  Bool_t MetaLHC16  = false;
  Bool_t MetaLHC17  = false;
  Bool_t MetaLHC18  = false;
  Bool_t LHC18q	    = false;
  Bool_t LHC18r	    = false;
  Bool_t LHC20g7a   = false;
  Bool_t LHC20g7b   = false;
  Bool_t LHC22f3    = false;

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


  Bool_t isProton	= false;
  Bool_t isDeuteron     = false;
  Bool_t isAntiProton   = false;
  Bool_t isAntiDeuteron = false;
  Bool_t isPion		= false;
  Bool_t isAntiPion     = false;

  if(ParticleSpecies == 1) isProton = true;
  if(ParticleSpecies == 2) isDeuteron = true;
  if(ParticleSpecies == 3) isAntiProton = true;
  if(ParticleSpecies == 4) isAntiDeuteron = true;
  if(ParticleSpecies == 5) isPion = true;
  if(ParticleSpecies == 6) isAntiPion = true;

  Double_t p = Track.P();
  Double_t Mass = 0.0;
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



  Double_t mean = Mean->Eval(p);
  Mean->Delete();

  Double_t Resolution = 0.0;

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

  Double_t ScaleFactor = 1.0-(Resolution);
  Double_t sigma = (mean*ScaleFactor) - mean;
  if(TMath::Abs(sigma) < 0.0001) return -999.0;

  SigmaParticle = (mean - SignalITS) / (sigma);
  if(TMath::IsNaN(SigmaParticle)) return -999.0;

  return SigmaParticle;

} // end of CalculateSigmadEdxITS








Float_t AliAnalysisTask_Ld_CreateTrees_PairsOnly::CalculateInvariantMassLambda(double Momentum1[3], double Momentum2[3], Int_t WhichMassHypothesis){

  Bool_t CalculateLambda      = false;
  Bool_t CalculateWrongLambda = false;
  Bool_t CalculateKaonShort   = false;
  Bool_t CalculatePhoton      = false;

  if(WhichMassHypothesis == 1) CalculateLambda	    = true;
  if(WhichMassHypothesis == 2) CalculateWrongLambda  = true;
  if(WhichMassHypothesis == 3) CalculateKaonShort    = true;
  if(WhichMassHypothesis == 4) CalculatePhoton	    = true;

  const Double_t MassProton   = 0.938272088;    // GeV/c²
  const Double_t MassPion     = 0.13957039;     // GeV/c²
  const Double_t MassElectron = 0.00051099895;  // GeV/c²

  Double_t m1 = -999.0;
  Double_t m2 = -999.0;
  Double_t p1x = -999.0;
  Double_t p1y = -999.0;
  Double_t p1z = -999.0;
  Double_t p2x = -999.0;
  Double_t p2y = -999.0;
  Double_t p2z = -999.0;


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

  if(CalculateWrongLambda == true){

    m1 = MassPion;
    p1x = Momentum1[0];
    p1y = Momentum1[1];
    p1z = Momentum1[2];

    m2	= MassProton;
    p2x = Momentum2[0];
    p2y = Momentum2[1];
    p2z = Momentum2[2];

  }

  if(CalculateKaonShort == true){

    m1 = MassPion;
    p1x = Momentum1[0];
    p1y = Momentum1[1];
    p1z = Momentum1[2];

    m2 = MassPion;
    p2x = Momentum2[0];
    p2y = Momentum2[1];
    p2z = Momentum2[2];

  }

  if(CalculatePhoton == true){

    m1 = MassElectron;
    p1x = Momentum1[0];
    p1y = Momentum1[1];
    p1z = Momentum1[2];

    m2 = MassElectron;
    p2x = Momentum2[0];
    p2y = Momentum2[1];
    p2z = Momentum2[2];

  }

  Double_t E1 = TMath::Sqrt((m1*m1) + (p1x*p1x) + (p1y*p1y) + (p1z*p1z));
  Double_t E2 = TMath::Sqrt((m2*m2) + (p2x*p2x) + (p2y*p2y) + (p2z*p2z));

  Float_t MassInvariant = TMath::Sqrt(((E1+E2)*(E1+E2)) - ((p1x+p2x)*(p1x+p2x)) - ((p1y+p2y)*(p1y+p2y)) - ((p1z+p2z)*(p1z+p2z)));

  return MassInvariant;

} // end of CalculateInvariantMassLambda












