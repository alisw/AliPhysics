#include "TChain.h" 
#include "TTree.h"
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
#include "TDatabasePDG.h"
#include "AliMultSelection.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODTrackSelection.h"
#include "AliVAODHeader.h"

#include "AliKFParticleBase.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <string>

#include "TDatabasePDG.h"
#include "TPDGCode.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTask_pd_CreateTrees_PairsOnly.h"

using namespace std;
ClassImp(AliAnalysisTask_pd_CreateTrees_PairsOnly) 






AliAnalysisTask_pd_CreateTrees_PairsOnly::AliAnalysisTask_pd_CreateTrees_PairsOnly() : AliAnalysisTaskSE(),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(0),
  fSaveTree_Proton(0),
  fProton_px(0),
  fProton_py(0),
  fProton_pz(0),
  fProton_pTPC(0),
  fProton_Eta(0),
  fProton_Phi(0),
  fProton_TPC_Chi2overNDF(0),
  fProton_TPC_dEdx(0),
  fProton_TPC_dEdx_Sigma(0),
  fProton_TOF_Beta(0),
  fProton_TOF_Beta_Sigma(0),
  fProton_TOF_Mass2(0),
  fProton_TOF_Mass2_Sigma(0),
  fProton_ITS_dEdx(0),
  fProton_ITS_dEdx_Sigma(0),
  fProton_DCAxy(0),
  fProton_DCAz(0),
  fProton_Event_Centrality(0),
  fProton_Event_PrimaryVertexZ(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nSharedCluster(0),
  fProton_TPC_nClusterFindable(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_nParticlesPerEvent(0),
  fProton_ID(0),
  fProton_Event_Identifier(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(0),
  fDeuteron_py(0),
  fDeuteron_pz(0),
  fDeuteron_pTPC(0),
  fDeuteron_Eta(0),
  fDeuteron_Phi(0),
  fDeuteron_TPC_Chi2overNDF(0),
  fDeuteron_TPC_dEdx(0),
  fDeuteron_TPC_dEdx_Sigma(0),
  fDeuteron_TOF_Beta(0),
  fDeuteron_TOF_Beta_Sigma(0),
  fDeuteron_TOF_Mass2(0),
  fDeuteron_TOF_Mass2_Sigma(0),
  fDeuteron_ITS_dEdx(0),
  fDeuteron_ITS_dEdx_Sigma(0),
  fDeuteron_DCAxy(0),
  fDeuteron_DCAz(0),
  fDeuteron_Event_Centrality(0),
  fDeuteron_Event_PrimaryVertexZ(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nSharedCluster(0),
  fDeuteron_TPC_nClusterFindable(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_nParticlesPerEvent(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Identifier(0),
  fSaveTree_AntiProton(0),
  fAntiProton_px(0),
  fAntiProton_py(0),
  fAntiProton_pz(0),
  fAntiProton_pTPC(0),
  fAntiProton_Eta(0),
  fAntiProton_Phi(0),
  fAntiProton_TPC_Chi2overNDF(0),
  fAntiProton_TPC_dEdx(0),
  fAntiProton_TPC_dEdx_Sigma(0),
  fAntiProton_TOF_Beta(0),
  fAntiProton_TOF_Beta_Sigma(0),
  fAntiProton_TOF_Mass2(0),
  fAntiProton_TOF_Mass2_Sigma(0),
  fAntiProton_ITS_dEdx(0),
  fAntiProton_ITS_dEdx_Sigma(0),
  fAntiProton_DCAxy(0),
  fAntiProton_DCAz(0),
  fAntiProton_Event_Centrality(0),
  fAntiProton_Event_PrimaryVertexZ(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nSharedCluster(0),
  fAntiProton_TPC_nClusterFindable(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_nParticlesPerEvent(0),
  fAntiProton_ID(0),
  fAntiProton_Event_Identifier(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(0),
  fAntiDeuteron_py(0),
  fAntiDeuteron_pz(0),
  fAntiDeuteron_pTPC(0),
  fAntiDeuteron_Eta(0),
  fAntiDeuteron_Phi(0),
  fAntiDeuteron_TPC_Chi2overNDF(0),
  fAntiDeuteron_TPC_dEdx(0),
  fAntiDeuteron_TPC_dEdx_Sigma(0),
  fAntiDeuteron_TOF_Beta(0),
  fAntiDeuteron_TOF_Beta_Sigma(0),
  fAntiDeuteron_TOF_Mass2(0),
  fAntiDeuteron_TOF_Mass2_Sigma(0),
  fAntiDeuteron_ITS_dEdx(0),
  fAntiDeuteron_ITS_dEdx_Sigma(0),
  fAntiDeuteron_DCAxy(0),
  fAntiDeuteron_DCAz(0),
  fAntiDeuteron_Event_Centrality(0),
  fAntiDeuteron_Event_PrimaryVertexZ(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nSharedCluster(0),
  fAntiDeuteron_TPC_nClusterFindable(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_nParticlesPerEvent(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Identifier(0)
{


}



AliAnalysisTask_pd_CreateTrees_PairsOnly::AliAnalysisTask_pd_CreateTrees_PairsOnly(const char *name,int CollisionSystem) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(CollisionSystem),
  fSaveTree_Proton(0),
  fProton_px(0),
  fProton_py(0),
  fProton_pz(0),
  fProton_pTPC(0),
  fProton_Eta(0),
  fProton_Phi(0),
  fProton_TPC_Chi2overNDF(0),
  fProton_TPC_dEdx(0),
  fProton_TPC_dEdx_Sigma(0),
  fProton_TOF_Beta(0),
  fProton_TOF_Beta_Sigma(0),
  fProton_TOF_Mass2(0),
  fProton_TOF_Mass2_Sigma(0),
  fProton_ITS_dEdx(0),
  fProton_ITS_dEdx_Sigma(0),
  fProton_DCAxy(0),
  fProton_DCAz(0),
  fProton_Event_Centrality(0),
  fProton_Event_PrimaryVertexZ(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nSharedCluster(0),
  fProton_TPC_nClusterFindable(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_nParticlesPerEvent(0),
  fProton_ID(0),
  fProton_Event_Identifier(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(0),
  fDeuteron_py(0),
  fDeuteron_pz(0),
  fDeuteron_pTPC(0),
  fDeuteron_Eta(0),
  fDeuteron_Phi(0),
  fDeuteron_TPC_Chi2overNDF(0),
  fDeuteron_TPC_dEdx(0),
  fDeuteron_TPC_dEdx_Sigma(0),
  fDeuteron_TOF_Beta(0),
  fDeuteron_TOF_Beta_Sigma(0),
  fDeuteron_TOF_Mass2(0),
  fDeuteron_TOF_Mass2_Sigma(0),
  fDeuteron_ITS_dEdx(0),
  fDeuteron_ITS_dEdx_Sigma(0),
  fDeuteron_DCAxy(0),
  fDeuteron_DCAz(0),
  fDeuteron_Event_Centrality(0),
  fDeuteron_Event_PrimaryVertexZ(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nSharedCluster(0),
  fDeuteron_TPC_nClusterFindable(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_nParticlesPerEvent(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Identifier(0),
  fSaveTree_AntiProton(0),
  fAntiProton_px(0),
  fAntiProton_py(0),
  fAntiProton_pz(0),
  fAntiProton_pTPC(0),
  fAntiProton_Eta(0),
  fAntiProton_Phi(0),
  fAntiProton_TPC_Chi2overNDF(0),
  fAntiProton_TPC_dEdx(0),
  fAntiProton_TPC_dEdx_Sigma(0),
  fAntiProton_TOF_Beta(0),
  fAntiProton_TOF_Beta_Sigma(0),
  fAntiProton_TOF_Mass2(0),
  fAntiProton_TOF_Mass2_Sigma(0),
  fAntiProton_ITS_dEdx(0),
  fAntiProton_ITS_dEdx_Sigma(0),
  fAntiProton_DCAxy(0),
  fAntiProton_DCAz(0),
  fAntiProton_Event_Centrality(0),
  fAntiProton_Event_PrimaryVertexZ(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nSharedCluster(0),
  fAntiProton_TPC_nClusterFindable(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_nParticlesPerEvent(0),
  fAntiProton_ID(0),
  fAntiProton_Event_Identifier(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(0),
  fAntiDeuteron_py(0),
  fAntiDeuteron_pz(0),
  fAntiDeuteron_pTPC(0),
  fAntiDeuteron_Eta(0),
  fAntiDeuteron_Phi(0),
  fAntiDeuteron_TPC_Chi2overNDF(0),
  fAntiDeuteron_TPC_dEdx(0),
  fAntiDeuteron_TPC_dEdx_Sigma(0),
  fAntiDeuteron_TOF_Beta(0),
  fAntiDeuteron_TOF_Beta_Sigma(0),
  fAntiDeuteron_TOF_Mass2(0),
  fAntiDeuteron_TOF_Mass2_Sigma(0),
  fAntiDeuteron_ITS_dEdx(0),
  fAntiDeuteron_ITS_dEdx_Sigma(0),
  fAntiDeuteron_DCAxy(0),
  fAntiDeuteron_DCAz(0),
  fAntiDeuteron_Event_Centrality(0),
  fAntiDeuteron_Event_PrimaryVertexZ(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nSharedCluster(0),
  fAntiDeuteron_TPC_nClusterFindable(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_nParticlesPerEvent(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Identifier(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());

}

  
AliAnalysisTask_pd_CreateTrees_PairsOnly::~AliAnalysisTask_pd_CreateTrees_PairsOnly()
{

  if(fSaveTree_Proton)
    {
      delete fSaveTree_Proton;
    }

  if(fSaveTree_Deuteron)
    {
      delete fSaveTree_Deuteron;
    }

  if(fSaveTree_AntiProton)
    {
      delete fSaveTree_AntiProton;
    }

  if(fSaveTree_AntiDeuteron)
    {
      delete fSaveTree_AntiDeuteron;
    }


}







void AliAnalysisTask_pd_CreateTrees_PairsOnly::UserCreateOutputObjects()
{

  fSaveTree_Proton = new TTree("fSaveTree_Proton","fSaveTree_Proton");
  fSaveTree_Proton->Branch("Proton_px",&fProton_px,"Proton_px/f");
  fSaveTree_Proton->Branch("Proton_py",&fProton_py,"Proton_py/f");
  fSaveTree_Proton->Branch("Proton_pz",&fProton_pz,"Proton_pz/f");
  fSaveTree_Proton->Branch("Proton_pTPC",&fProton_pTPC,"Proton_pTPC/f");
  fSaveTree_Proton->Branch("Proton_Eta",&fProton_Eta,"Proton_Eta/f");
  fSaveTree_Proton->Branch("Proton_Phi",&fProton_Phi,"Proton_Phi/f");
  fSaveTree_Proton->Branch("Proton_TPC_Chi2overNDF",&fProton_TPC_Chi2overNDF,"Proton_TPC_Chi2overNDF/f");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx",&fProton_TPC_dEdx,"Proton_TPC_dEdx/f");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx_Sigma",&fProton_TPC_dEdx_Sigma,"Proton_TPC_dEdx_Sigma/f");
  fSaveTree_Proton->Branch("Proton_TOF_Beta",&fProton_TOF_Beta,"Proton_TOF_Beta/f");
  fSaveTree_Proton->Branch("Proton_TOF_Beta_Sigma",&fProton_TOF_Beta_Sigma,"Proton_TOF_Beta_Sigma/f");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2",&fProton_TOF_Mass2,"Proton_TOF_Mass2/f");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2_Sigma",&fProton_TOF_Mass2_Sigma,"Proton_TOF_Mass2_Sigma/f");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx",&fProton_ITS_dEdx,"Proton_ITS_dEdx/f");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx_Sigma",&fProton_ITS_dEdx_Sigma,"Proton_ITS_dEdx_Sigma/f");
  fSaveTree_Proton->Branch("Proton_DCAxy",&fProton_DCAxy,"Proton_DCAxy/f");
  fSaveTree_Proton->Branch("Proton_DCAz",&fProton_DCAz,"Proton_DCAz/f");
  fSaveTree_Proton->Branch("Proton_Event_Centrality",&fProton_Event_Centrality,"Proton_Event_Centrality/f");
  fSaveTree_Proton->Branch("Proton_Event_PrimaryVertexZ",&fProton_Event_PrimaryVertexZ,"Proton_Event_PrimaryVertexZ/f");
  fSaveTree_Proton->Branch("Proton_TPC_nCrossedRows",&fProton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fSaveTree_Proton->Branch("Proton_TPC_nSharedCluster",&fProton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fSaveTree_Proton->Branch("Proton_TPC_nClusterFindable",&fProton_TPC_nClusterFindable,"Proton_TPC_nClusterFindable/s");
  fSaveTree_Proton->Branch("Proton_TPC_nCluster",&fProton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fSaveTree_Proton->Branch("Proton_ITS_nCluster",&fProton_ITS_nCluster,"Proton_ITS_nCluster/s");
  fSaveTree_Proton->Branch("Proton_nParticlesPerEvent",&fProton_nParticlesPerEvent,"Proton_nParticlesPerEvent/s");
  fSaveTree_Proton->Branch("Proton_ID",&fProton_ID,"Proton_ID/i");
  fSaveTree_Proton->Branch("Proton_Event_Identifier",&fProton_Event_Identifier,"Proton_Event_Identifier/i");


  fSaveTree_Deuteron = new TTree("fSaveTree_Deuteron","fSaveTree_Deuteron");
  fSaveTree_Deuteron->Branch("Deuteron_px",&fDeuteron_px,"Deuteron_px/f");
  fSaveTree_Deuteron->Branch("Deuteron_py",&fDeuteron_py,"Deuteron_py/f");
  fSaveTree_Deuteron->Branch("Deuteron_pz",&fDeuteron_pz,"Deuteron_pz/f");
  fSaveTree_Deuteron->Branch("Deuteron_pTPC",&fDeuteron_pTPC,"Deuteron_pTPC/f");
  fSaveTree_Deuteron->Branch("Deuteron_Eta",&fDeuteron_Eta,"Deuteron_Eta/f");
  fSaveTree_Deuteron->Branch("Deuteron_Phi",&fDeuteron_Phi,"Deuteron_Phi/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_Chi2overNDF",&fDeuteron_TPC_Chi2overNDF,"Deuteron_TPC_Chi2overNDF/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx",&fDeuteron_TPC_dEdx,"Deuteron_TPC_dEdx/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx_Sigma",&fDeuteron_TPC_dEdx_Sigma,"Deuteron_TPC_dEdx_Sigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Beta",&fDeuteron_TOF_Beta,"Deuteron_TOF_Beta/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Beta_Sigma",&fDeuteron_TOF_Beta_Sigma,"Deuteron_TOF_Beta_Sigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2",&fDeuteron_TOF_Mass2,"Deuteron_TOF_Mass2/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2_Sigma",&fDeuteron_TOF_Mass2_Sigma,"Deuteron_TOF_Mass2_Sigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx",&fDeuteron_ITS_dEdx,"Deuteron_ITS_dEdx/f");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx_Sigma",&fDeuteron_ITS_dEdx_Sigma,"Deuteron_ITS_dEdx_Sigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_DCAxy",&fDeuteron_DCAxy,"Deuteron_DCAxy/f");
  fSaveTree_Deuteron->Branch("Deuteron_DCAz",&fDeuteron_DCAz,"Deuteron_DCAz/f");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Centrality",&fDeuteron_Event_Centrality,"Deuteron_Event_Centrality/f");
  fSaveTree_Deuteron->Branch("Deuteron_Event_PrimaryVertexZ",&fDeuteron_Event_PrimaryVertexZ,"Deuteron_Event_PrimaryVertexZ/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&fDeuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&fDeuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nClusterFindable",&fDeuteron_TPC_nClusterFindable,"Deuteron_TPC_nClusterFindable/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCluster",&fDeuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_nCluster",&fDeuteron_ITS_nCluster,"Deuteron_ITS_nCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_nParticlesPerEvent",&fDeuteron_nParticlesPerEvent,"Deuteron_nParticlesPerEvent/s");
  fSaveTree_Deuteron->Branch("Deuteron_ID",&fDeuteron_ID,"Deuteron_ID/i");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Identifier",&fDeuteron_Event_Identifier,"Deuteron_Event_Identifier/i");




  fSaveTree_AntiProton = new TTree("fSaveTree_AntiProton","fSaveTree_AntiProton");
  fSaveTree_AntiProton->Branch("AntiProton_px",&fAntiProton_px,"AntiProton_px/f");
  fSaveTree_AntiProton->Branch("AntiProton_py",&fAntiProton_py,"AntiProton_py/f");
  fSaveTree_AntiProton->Branch("AntiProton_pz",&fAntiProton_pz,"AntiProton_pz/f");
  fSaveTree_AntiProton->Branch("AntiProton_pTPC",&fAntiProton_pTPC,"AntiProton_pTPC/f");
  fSaveTree_AntiProton->Branch("AntiProton_Eta",&fAntiProton_Eta,"AntiProton_Eta/f");
  fSaveTree_AntiProton->Branch("AntiProton_Phi",&fAntiProton_Phi,"AntiProton_Phi/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_Chi2overNDF",&fAntiProton_TPC_Chi2overNDF,"AntiProton_TPC_Chi2overNDF/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx",&fAntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx_Sigma",&fAntiProton_TPC_dEdx_Sigma,"AntiProton_TPC_dEdx_Sigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Beta",&fAntiProton_TOF_Beta,"AntiProton_TOF_Beta/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Beta_Sigma",&fAntiProton_TOF_Beta_Sigma,"AntiProton_TOF_Beta_Sigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2",&fAntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2_Sigma",&fAntiProton_TOF_Mass2_Sigma,"AntiProton_TOF_Mass2_Sigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx",&fAntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/f");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx_Sigma",&fAntiProton_ITS_dEdx_Sigma,"AntiProton_ITS_dEdx_Sigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_DCAxy",&fAntiProton_DCAxy,"AntiProton_DCAxy/f");
  fSaveTree_AntiProton->Branch("AntiProton_DCAz",&fAntiProton_DCAz,"AntiProton_DCAz/f");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Centrality",&fAntiProton_Event_Centrality,"AntiProton_Event_Centrality/f");
  fSaveTree_AntiProton->Branch("AntiProton_Event_PrimaryVertexZ",&fAntiProton_Event_PrimaryVertexZ,"AntiProton_Event_PrimaryVertexZ/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&fAntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&fAntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nClusterFindable",&fAntiProton_TPC_nClusterFindable,"AntiProton_TPC_nClusterFindable/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCluster",&fAntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_nCluster",&fAntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_nParticlesPerEvent",&fAntiProton_nParticlesPerEvent,"AntiProton_nParticlesPerEvent/s");
  fSaveTree_AntiProton->Branch("AntiProton_ID",&fAntiProton_ID,"AntiProton_ID/i");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Identifier",&fAntiProton_Event_Identifier,"AntiProton_Event_Identifier/i");


  fSaveTree_AntiDeuteron = new TTree("fSaveTree_AntiDeuteron","fSaveTree_AntiDeuteron");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_px",&fAntiDeuteron_px,"AntiDeuteron_px/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_py",&fAntiDeuteron_py,"AntiDeuteron_py/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pz",&fAntiDeuteron_pz,"AntiDeuteron_pz/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&fAntiDeuteron_pTPC,"AntiDeuteron_pTPC/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&fAntiDeuteron_Eta,"AntiDeuteron_Eta/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&fAntiDeuteron_Phi,"AntiDeuteron_Phi/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2overNDF",&fAntiDeuteron_TPC_Chi2overNDF,"AntiDeuteron_TPC_Chi2overNDF/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&fAntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_Sigma",&fAntiDeuteron_TPC_dEdx_Sigma,"AntiDeuteron_TPC_dEdx_Sigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta",&fAntiDeuteron_TOF_Beta,"AntiDeuteron_TOF_Beta/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta_Sigma",&fAntiDeuteron_TOF_Beta_Sigma,"AntiDeuteron_TOF_Beta_Sigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&fAntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_Sigma",&fAntiDeuteron_TOF_Mass2_Sigma,"AntiDeuteron_TOF_Mass2_Sigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&fAntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_Sigma",&fAntiDeuteron_ITS_dEdx_Sigma,"AntiDeuteron_ITS_dEdx_Sigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&fAntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&fAntiDeuteron_DCAz,"AntiDeuteron_DCAz/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Centrality",&fAntiDeuteron_Event_Centrality,"AntiDeuteron_Event_Centrality/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_PrimaryVertexZ",&fAntiDeuteron_Event_PrimaryVertexZ,"AntiDeuteron_Event_PrimaryVertexZ/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&fAntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&fAntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nClusterFindable",&fAntiDeuteron_TPC_nClusterFindable,"AntiDeuteron_TPC_nClusterFindable/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&fAntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&fAntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_nParticlesPerEvent",&fAntiDeuteron_nParticlesPerEvent,"AntiDeuteron_nParticlesPerEvent/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ID",&fAntiDeuteron_ID,"AntiDeuteron_ID/i");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&fAntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/i");




  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiProton);
  PostData(4,fSaveTree_AntiDeuteron);



} // end of UserCreateOutputObjects









void AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec(Option_t*)
{

  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAODEvent)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AOD event found!");

  fHeader = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  if(!fHeader)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No Header found!");

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No PIDResponse found!");


  // debug analysis
  bool DebugEventSelection  = false;


 
  bool isPbPb = true;



  // define event cuts
  double PrimaryVertexMaxZ  = 10.0; // cm
  double Centrality_min	    = 0.0;
  double Centrality_max	    = 100.0;

  if(fCollisionSystem == 1)
  {
    Centrality_min = 0.0;
    Centrality_max = 10.0;
  }

  if(fCollisionSystem == 2)
  {
    Centrality_min = 30.0;
    Centrality_max = 50.0;
  }


  // use only events containing tracks
  int nTracks = fAODEvent->GetNumberOfTracks();
  if(nTracks == 0) return;
  
  // get primary vertex
  AliAODVertex *PrimaryVertex = fAODEvent->GetPrimaryVertex();
  if(!PrimaryVertex)::Warning("AliAnalsisTask_pd_CreateTrees_PairsOnlyd::UserExec","No AliAODVertex object found!");
  double PrimaryVertexPos[3] = {-999.0,-999.0,-999.0};
  PrimaryVertex->GetXYZ(PrimaryVertexPos);

  // apply cut on z-position of primary vertex
  double PrimaryVertexZ = PrimaryVertexPos[2]; // cm
  if(TMath::Abs(PrimaryVertexZ) > PrimaryVertexMaxZ) return;

  // apply centrality cut
  double Centrality = -999.0;
  if(isPbPb == true){

    AliMultSelection *MultSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
    Centrality = MultSelection->GetMultiplicityPercentile("V0M");
    if((Centrality < Centrality_min) || (Centrality > Centrality_max)) return;

  }




  // get event information
  unsigned int PeriodNumber	= fAODEvent->GetPeriodNumber();
  unsigned int OrbitNumber	= fAODEvent->GetOrbitNumber();
  unsigned int BunchCrossNumber	= fAODEvent->GetBunchCrossNumber();

  // print event information
  if(DebugEventSelection)
  {

    cout << "" << endl;
    cout << "PeriodNumber:\t\t\t" << PeriodNumber << endl;
    cout << "OrbitNumber:\t\t\t" << OrbitNumber << endl;
    cout << "BunchCrossNumber:\t\t" << BunchCrossNumber << endl;
    cout << "isPbPb:\t\t\t\t" << isPbPb << endl;
    cout << "Centrality:\t\t\t" << Centrality << " %" << endl;
    cout << "z-position of primary vertex:\t" << PrimaryVertexZ << " cm" << endl;
    cout << "Number of tracks in event:\t" << nTracks << endl;

  } // end of DebugEventSelection





  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ proton selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  float     Proton_px;
  float     Proton_py;
  float     Proton_pz;
  float     Proton_pTPC;
  float     Proton_Eta;
  float     Proton_Phi;
  float     Proton_TPC_Chi2overNDF;
  float     Proton_TPC_dEdx;
  float     Proton_TPC_dEdx_Sigma;
  float     Proton_TOF_Beta;
  float     Proton_TOF_Beta_Sigma;
  float     Proton_TOF_Mass2;
  float     Proton_TOF_Mass2_Sigma;
  float     Proton_ITS_dEdx;
  float     Proton_ITS_dEdx_Sigma;
  float     Proton_DCAxy;
  float     Proton_DCAz;
  float     Proton_Event_Centrality;
  float     Proton_Event_PrimaryVertexZ;
  unsigned short    Proton_TPC_nCrossedRows;
  unsigned short    Proton_TPC_nSharedCluster;
  unsigned short    Proton_TPC_nClusterFindable;
  unsigned short    Proton_TPC_nCluster;
  unsigned short    Proton_ITS_nCluster;
  unsigned int      Proton_ID;
  unsigned int      Proton_Event_Identifier;

  TTree *fTempTree_Proton = new TTree("fTempTree_Proton","fTempTree_Proton");
  fTempTree_Proton->Branch("Proton_px",&Proton_px,"Proton_px/f");
  fTempTree_Proton->Branch("Proton_py",&Proton_py,"Proton_py/f");
  fTempTree_Proton->Branch("Proton_pz",&Proton_pz,"Proton_pz/f");
  fTempTree_Proton->Branch("Proton_pTPC",&Proton_pTPC,"Proton_pTPC/f");
  fTempTree_Proton->Branch("Proton_Eta",&Proton_Eta,"Proton_Eta/f");
  fTempTree_Proton->Branch("Proton_Phi",&Proton_Phi,"Proton_Phi/f");
  fTempTree_Proton->Branch("Proton_TPC_Chi2overNDF",&Proton_TPC_Chi2overNDF,"Proton_TPC_Chi2overNDF/f");
  fTempTree_Proton->Branch("Proton_TPC_dEdx",&Proton_TPC_dEdx,"Proton_TPC_dEdx/f");
  fTempTree_Proton->Branch("Proton_TPC_dEdx_Sigma",&Proton_TPC_dEdx_Sigma,"Proton_TPC_dEdx_Sigma/f");
  fTempTree_Proton->Branch("Proton_TOF_Beta",&Proton_TOF_Beta,"Proton_TOF_Beta/f");
  fTempTree_Proton->Branch("Proton_TOF_Beta_Sigma",&Proton_TOF_Beta_Sigma,"Proton_TOF_Beta_Sigma/f");
  fTempTree_Proton->Branch("Proton_TOF_Mass2",&Proton_TOF_Mass2,"Proton_TOF_Mass2/f");
  fTempTree_Proton->Branch("Proton_TOF_Mass2_Sigma",&Proton_TOF_Mass2_Sigma,"Proton_TOF_Mass2_Sigma/f");
  fTempTree_Proton->Branch("Proton_ITS_dEdx",&Proton_ITS_dEdx,"Proton_ITS_dEdx/f");
  fTempTree_Proton->Branch("Proton_ITS_dEdx_Sigma",&Proton_ITS_dEdx_Sigma,"Proton_ITS_dEdx_Sigma/f");
  fTempTree_Proton->Branch("Proton_DCAxy",&Proton_DCAxy,"Proton_DCAxy/f");
  fTempTree_Proton->Branch("Proton_DCAz",&Proton_DCAz,"Proton_DCAz/f");
  fTempTree_Proton->Branch("Proton_Event_Centrality",&Proton_Event_Centrality,"Proton_Event_Centrality/f");
  fTempTree_Proton->Branch("Proton_Event_PrimaryVertexZ",&Proton_Event_PrimaryVertexZ,"Proton_Event_PrimaryVertexZ/f");
  fTempTree_Proton->Branch("Proton_TPC_nCrossedRows",&Proton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fTempTree_Proton->Branch("Proton_TPC_nSharedCluster",&Proton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fTempTree_Proton->Branch("Proton_TPC_nClusterFindable",&Proton_TPC_nClusterFindable,"Proton_TPC_nClusterFindable/s");
  fTempTree_Proton->Branch("Proton_TPC_nCluster",&Proton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fTempTree_Proton->Branch("Proton_ITS_nCluster",&Proton_ITS_nCluster,"Proton_ITS_nCluster/s");
  fTempTree_Proton->Branch("Proton_ID",&Proton_ID,"Proton_ID/i");
  fTempTree_Proton->Branch("Proton_Event_Identifier",&Proton_Event_Identifier,"Proton_Event_Identifier/i");

  unsigned short nProtonsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// proton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply proton cuts
    bool PassedProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,true);
    if(!PassedProtonCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    Proton_px			    = Track->Px();
    Proton_py			    = Track->Py();
    Proton_pz			    = Track->Pz();
    Proton_pTPC			    = Track->GetTPCmomentum();
    Proton_Eta			    = Track->Eta();
    Proton_Phi			    = Track->Phi();
    Proton_TPC_Chi2overNDF	    = Track->GetTPCchi2perNDF();
    Proton_TPC_dEdx		    = Track->GetTPCsignal();
    Proton_TPC_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    Proton_TOF_Beta		    = Track->GetTOFsignal();
    Proton_TOF_Beta_Sigma	    = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
    Proton_TOF_Mass2		    = CalculateMassSquareTOF(*Track);
    Proton_ITS_dEdx		    = Track->GetITSsignal();
    Proton_ITS_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kProton);
    Proton_DCAxy		    = DCAxy;
    Proton_DCAz			    = DCAz;
    Proton_Event_Centrality	    = Centrality;
    Proton_Event_PrimaryVertexZ    = PrimaryVertexZ;
    Proton_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    Proton_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    Proton_TPC_nClusterFindable	    = Track->GetTPCNclsF();
    Proton_TPC_nCluster		    = Track->GetTPCNcls();
    Proton_ITS_nCluster		    = Track->GetITSNcls();
    Proton_ID			    = Track->GetID();
    Proton_Event_Identifier	    = BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    fTempTree_Proton->Fill();
    nProtonsSelected++;

  } // end of proton loop




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ deuteron selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float     Deuteron_px;
  float     Deuteron_py;
  float     Deuteron_pz;
  float     Deuteron_pTPC;
  float     Deuteron_Eta;
  float     Deuteron_Phi;
  float     Deuteron_TPC_Chi2overNDF;
  float     Deuteron_TPC_dEdx;
  float     Deuteron_TPC_dEdx_Sigma;
  float     Deuteron_TOF_Beta;
  float     Deuteron_TOF_Beta_Sigma;
  float     Deuteron_TOF_Mass2;
  float     Deuteron_TOF_Mass2_Sigma;
  float     Deuteron_ITS_dEdx;
  float     Deuteron_ITS_dEdx_Sigma;
  float     Deuteron_DCAxy;
  float     Deuteron_DCAz;
  float     Deuteron_Event_Centrality;
  float     Deuteron_Event_PrimaryVertexZ;
  unsigned short  Deuteron_TPC_nCrossedRows;
  unsigned short  Deuteron_TPC_nSharedCluster;
  unsigned short  Deuteron_TPC_nClusterFindable;
  unsigned short  Deuteron_TPC_nCluster;
  unsigned short  Deuteron_ITS_nCluster;
  unsigned int    Deuteron_ID;
  unsigned int    Deuteron_Event_Identifier;

  TTree *fTempTree_Deuteron = new TTree("fTempTree_Deuteron","fTempTree_Deuteron");
  fTempTree_Deuteron->Branch("Deuteron_px",&Deuteron_px,"Deuteron_px/f");
  fTempTree_Deuteron->Branch("Deuteron_py",&Deuteron_py,"Deuteron_py/f");
  fTempTree_Deuteron->Branch("Deuteron_pz",&Deuteron_pz,"Deuteron_pz/f");
  fTempTree_Deuteron->Branch("Deuteron_pTPC",&Deuteron_pTPC,"Deuteron_pTPC/f");
  fTempTree_Deuteron->Branch("Deuteron_Eta",&Deuteron_Eta,"Deuteron_Eta/f");
  fTempTree_Deuteron->Branch("Deuteron_Phi",&Deuteron_Phi,"Deuteron_Phi/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_Chi2overNDF",&Deuteron_TPC_Chi2overNDF,"Deuteron_TPC_Chi2overNDF/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_dEdx",&Deuteron_TPC_dEdx,"Deuteron_TPC_dEdx/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_dEdx_Sigma",&Deuteron_TPC_dEdx_Sigma,"Deuteron_TPC_dEdx_Sigma/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Beta",&Deuteron_TOF_Beta,"Deuteron_TOF_Beta/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Beta_Sigma",&Deuteron_TOF_Beta_Sigma,"Deuteron_TOF_Beta_Sigma/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Mass2",&Deuteron_TOF_Mass2,"Deuteron_TOF_Mass2/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Mass2_Sigma",&Deuteron_TOF_Mass2_Sigma,"Deuteron_TOF_Mass2_Sigma/f");
  fTempTree_Deuteron->Branch("Deuteron_ITS_dEdx",&Deuteron_ITS_dEdx,"Deuteron_ITS_dEdx/f");
  fTempTree_Deuteron->Branch("Deuteron_ITS_dEdx_Sigma",&Deuteron_ITS_dEdx_Sigma,"Deuteron_ITS_dEdx_Sigma/f");
  fTempTree_Deuteron->Branch("Deuteron_DCAxy",&Deuteron_DCAxy,"Deuteron_DCAxy/f");
  fTempTree_Deuteron->Branch("Deuteron_DCAz",&Deuteron_DCAz,"Deuteron_DCAz/f");
  fTempTree_Deuteron->Branch("Deuteron_Event_Centrality",&Deuteron_Event_Centrality,"Deuteron_Event_Centrality/f");
  fTempTree_Deuteron->Branch("Deuteron_Event_PrimaryVertexZ",&Deuteron_Event_PrimaryVertexZ,"Deuteron_Event_PrimaryVertexZ/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&Deuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&Deuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nClusterFindable",&Deuteron_TPC_nClusterFindable,"Deuteron_TPC_nClusterFindable/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nCluster",&Deuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_ITS_nCluster",&Deuteron_ITS_nCluster,"Deuteron_ITS_nCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_ID",&Deuteron_ID,"Deuteron_ID/i");
  fTempTree_Deuteron->Branch("Deuteron_Event_Identifier",&Deuteron_Event_Identifier,"Deuteron_Event_Identifier/i");

  unsigned short nDeuteronsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// deuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply deuteron cuts
    bool PassedDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,true);
    if(!PassedDeuteronCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    Deuteron_px			    = Track->Px();
    Deuteron_py			    = Track->Py();
    Deuteron_pz			    = Track->Pz();
    Deuteron_pTPC		    = Track->GetTPCmomentum();
    Deuteron_Eta		    = Track->Eta();
    Deuteron_Phi		    = Track->Phi();
    Deuteron_TPC_Chi2overNDF	    = Track->GetTPCchi2perNDF();
    Deuteron_TPC_dEdx		    = Track->GetTPCsignal();
    Deuteron_TPC_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    Deuteron_TOF_Beta		    = Track->GetTOFsignal();
    Deuteron_TOF_Beta_Sigma	    = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
    Deuteron_TOF_Mass2		    = CalculateMassSquareTOF(*Track);
    Deuteron_TOF_Mass2_Sigma	    = CalculateDeuteronSigmaMassSquareTOF(Track->Pt(),Deuteron_TOF_Mass2,true);
    Deuteron_ITS_dEdx		    = Track->GetITSsignal();
    Deuteron_ITS_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kDeuteron);
    Deuteron_DCAxy		    = DCAxy;
    Deuteron_DCAz		    = DCAz;
    Deuteron_Event_Centrality	    = Centrality;
    Deuteron_Event_PrimaryVertexZ  = PrimaryVertexZ;
    Deuteron_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    Deuteron_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    Deuteron_TPC_nClusterFindable   = Track->GetTPCNclsF();
    Deuteron_TPC_nCluster	    = Track->GetTPCNcls();
    Deuteron_ITS_nCluster	    = Track->GetITSNcls();
    Deuteron_ID			    = Track->GetID();
    Deuteron_Event_Identifier	    = BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);

    fTempTree_Deuteron->Fill();
    nDeuteronsSelected++;

  } // end of deuteron loop



  if((nProtonsSelected > 0) && (nDeuteronsSelected > 0)){

    for(int Proton = 0; Proton < nProtonsSelected; Proton++){

      TBranch *Branch_Proton_px			  = fTempTree_Proton->GetBranch("Proton_px");
      TBranch *Branch_Proton_py			  = fTempTree_Proton->GetBranch("Proton_py");
      TBranch *Branch_Proton_pz			  = fTempTree_Proton->GetBranch("Proton_pz");
      TBranch *Branch_Proton_pTPC		  = fTempTree_Proton->GetBranch("Proton_pTPC");
      TBranch *Branch_Proton_Eta		  = fTempTree_Proton->GetBranch("Proton_Eta");
      TBranch *Branch_Proton_Phi		  = fTempTree_Proton->GetBranch("Proton_Phi");
      TBranch *Branch_Proton_TPC_Chi2overNDF	  = fTempTree_Proton->GetBranch("Proton_TPC_Chi2overNDF");
      TBranch *Branch_Proton_TPC_dEdx		  = fTempTree_Proton->GetBranch("Proton_TPC_dEdx");
      TBranch *Branch_Proton_TPC_dEdx_Sigma	  = fTempTree_Proton->GetBranch("Proton_TPC_dEdx_Sigma");
      TBranch *Branch_Proton_TOF_Beta		  = fTempTree_Proton->GetBranch("Proton_TOF_Beta");
      TBranch *Branch_Proton_TOF_Beta_Sigma	  = fTempTree_Proton->GetBranch("Proton_TOF_Beta_Sigma");
      TBranch *Branch_Proton_TOF_Mass2		  = fTempTree_Proton->GetBranch("Proton_TOF_Mass2");
      TBranch *Branch_Proton_TOF_Mass2_Sigma	  = fTempTree_Proton->GetBranch("Proton_TOF_Mass2_Sigma");
      TBranch *Branch_Proton_ITS_dEdx		  = fTempTree_Proton->GetBranch("Proton_ITS_dEdx");
      TBranch *Branch_Proton_ITS_dEdx_Sigma	  = fTempTree_Proton->GetBranch("Proton_ITS_dEdx_Sigma");
      TBranch *Branch_Proton_DCAxy		  = fTempTree_Proton->GetBranch("Proton_DCAxy");
      TBranch *Branch_Proton_DCAz		  = fTempTree_Proton->GetBranch("Proton_DCAz");
      TBranch *Branch_Proton_Event_Centrality	  = fTempTree_Proton->GetBranch("Proton_Event_Centrality");
      TBranch *Branch_Proton_Event_PrimaryVertexZ = fTempTree_Proton->GetBranch("Proton_Event_PrimaryVertexZ");
      TBranch *Branch_Proton_TPC_nCrossedRows	  = fTempTree_Proton->GetBranch("Proton_TPC_nCrossedRows");
      TBranch *Branch_Proton_TPC_nSharedCluster	  = fTempTree_Proton->GetBranch("Proton_TPC_nSharedCluster");
      TBranch *Branch_Proton_TPC_nClusterFindable = fTempTree_Proton->GetBranch("Proton_TPC_nClusterFindable");
      TBranch *Branch_Proton_TPC_nCluster	  = fTempTree_Proton->GetBranch("Proton_TPC_nCluster");
      TBranch *Branch_Proton_ITS_nCluster	  = fTempTree_Proton->GetBranch("Proton_ITS_nCluster");
      TBranch *Branch_Proton_ID			  = fTempTree_Proton->GetBranch("Proton_ID");
      TBranch *Branch_Proton_Event_Identifier	  = fTempTree_Proton->GetBranch("Proton_Event_Identifier");


      Branch_Proton_px->SetAddress(&fProton_px);
      Branch_Proton_py->SetAddress(&fProton_py);
      Branch_Proton_pz->SetAddress(&fProton_pz);
      Branch_Proton_pTPC->SetAddress(&fProton_pTPC);
      Branch_Proton_Eta->SetAddress(&fProton_Eta);
      Branch_Proton_Phi->SetAddress(&fProton_Phi);
      Branch_Proton_TPC_Chi2overNDF->SetAddress(&fProton_TPC_Chi2overNDF);
      Branch_Proton_TPC_dEdx->SetAddress(&fProton_TPC_dEdx);
      Branch_Proton_TPC_dEdx_Sigma->SetAddress(&fProton_TPC_dEdx_Sigma);
      Branch_Proton_TOF_Beta->SetAddress(&fProton_TOF_Beta);
      Branch_Proton_TOF_Beta_Sigma->SetAddress(&fProton_TOF_Beta_Sigma);
      Branch_Proton_TOF_Mass2->SetAddress(&fProton_TOF_Mass2);
      Branch_Proton_TOF_Mass2_Sigma->SetAddress(&fProton_TOF_Mass2_Sigma);
      Branch_Proton_ITS_dEdx->SetAddress(&fProton_ITS_dEdx);
      Branch_Proton_ITS_dEdx_Sigma->SetAddress(&fProton_ITS_dEdx_Sigma);
      Branch_Proton_DCAxy->SetAddress(&fProton_DCAxy);
      Branch_Proton_DCAz->SetAddress(&fProton_DCAz);
      Branch_Proton_Event_Centrality->SetAddress(&fProton_Event_Centrality);
      Branch_Proton_Event_PrimaryVertexZ->SetAddress(&fProton_Event_PrimaryVertexZ);
      Branch_Proton_TPC_nCrossedRows->SetAddress(&fProton_TPC_nCrossedRows);
      Branch_Proton_TPC_nSharedCluster->SetAddress(&fProton_TPC_nSharedCluster);
      Branch_Proton_TPC_nClusterFindable->SetAddress(&fProton_TPC_nClusterFindable);
      Branch_Proton_TPC_nCluster->SetAddress(&fProton_TPC_nCluster);
      Branch_Proton_ITS_nCluster->SetAddress(&fProton_ITS_nCluster);
      Branch_Proton_ID->SetAddress(&fProton_ID);
      Branch_Proton_Event_Identifier->SetAddress(&fProton_Event_Identifier);


      Branch_Proton_px->SetAutoDelete(true);
      Branch_Proton_py->SetAutoDelete(true);
      Branch_Proton_pz->SetAutoDelete(true);
      Branch_Proton_pTPC->SetAutoDelete(true);
      Branch_Proton_Eta->SetAutoDelete(true);
      Branch_Proton_Phi->SetAutoDelete(true);
      Branch_Proton_TPC_Chi2overNDF->SetAutoDelete(true);
      Branch_Proton_TPC_dEdx->SetAutoDelete(true);
      Branch_Proton_TPC_dEdx_Sigma->SetAutoDelete(true);
      Branch_Proton_TOF_Beta->SetAutoDelete(true);
      Branch_Proton_TOF_Beta_Sigma->SetAutoDelete(true);
      Branch_Proton_TOF_Mass2->SetAutoDelete(true);
      Branch_Proton_TOF_Mass2_Sigma->SetAutoDelete(true);
      Branch_Proton_ITS_dEdx->SetAutoDelete(true);
      Branch_Proton_ITS_dEdx_Sigma->SetAutoDelete(true);
      Branch_Proton_DCAxy->SetAutoDelete(true);
      Branch_Proton_DCAz->SetAutoDelete(true);
      Branch_Proton_Event_Centrality->SetAutoDelete(true);
      Branch_Proton_Event_PrimaryVertexZ->SetAutoDelete(true);
      Branch_Proton_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_Proton_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_Proton_TPC_nClusterFindable->SetAutoDelete(true);
      Branch_Proton_TPC_nCluster->SetAutoDelete(true);
      Branch_Proton_ITS_nCluster->SetAutoDelete(true);
      Branch_Proton_ID->SetAutoDelete(true);
      Branch_Proton_Event_Identifier->SetAutoDelete(true);

      Branch_Proton_px->GetEntry(Proton);
      Branch_Proton_py->GetEntry(Proton);
      Branch_Proton_pz->GetEntry(Proton);
      Branch_Proton_pTPC->GetEntry(Proton);
      Branch_Proton_Eta->GetEntry(Proton);
      Branch_Proton_Phi->GetEntry(Proton);
      Branch_Proton_TPC_Chi2overNDF->GetEntry(Proton);
      Branch_Proton_TPC_dEdx->GetEntry(Proton);
      Branch_Proton_TPC_dEdx_Sigma->GetEntry(Proton);
      Branch_Proton_TOF_Beta->GetEntry(Proton);
      Branch_Proton_TOF_Beta_Sigma->GetEntry(Proton);
      Branch_Proton_TOF_Mass2->GetEntry(Proton);
      Branch_Proton_TOF_Mass2_Sigma->GetEntry(Proton);
      Branch_Proton_ITS_dEdx->GetEntry(Proton);
      Branch_Proton_ITS_dEdx_Sigma->GetEntry(Proton);
      Branch_Proton_DCAxy->GetEntry(Proton);
      Branch_Proton_DCAz->GetEntry(Proton);
      Branch_Proton_Event_Centrality->GetEntry(Proton);
      Branch_Proton_Event_PrimaryVertexZ->GetEntry(Proton);
      Branch_Proton_TPC_nCrossedRows->GetEntry(Proton);
      Branch_Proton_TPC_nSharedCluster->GetEntry(Proton);
      Branch_Proton_TPC_nClusterFindable->GetEntry(Proton);
      Branch_Proton_TPC_nCluster->GetEntry(Proton);
      Branch_Proton_ITS_nCluster->GetEntry(Proton);
      Branch_Proton_ID->GetEntry(Proton);
      Branch_Proton_Event_Identifier->GetEntry(Proton);
      fProton_nParticlesPerEvent = nProtonsSelected;

      fSaveTree_Proton->Fill();

    } // end of loop (copy protons)


    for(int Deuteron = 0; Deuteron < nDeuteronsSelected; Deuteron++){

      TBranch *Branch_Deuteron_px		    = fTempTree_Deuteron->GetBranch("Deuteron_px");
      TBranch *Branch_Deuteron_py		    = fTempTree_Deuteron->GetBranch("Deuteron_py");
      TBranch *Branch_Deuteron_pz		    = fTempTree_Deuteron->GetBranch("Deuteron_pz");
      TBranch *Branch_Deuteron_pTPC		    = fTempTree_Deuteron->GetBranch("Deuteron_pTPC");
      TBranch *Branch_Deuteron_Eta		    = fTempTree_Deuteron->GetBranch("Deuteron_Eta");
      TBranch *Branch_Deuteron_Phi		    = fTempTree_Deuteron->GetBranch("Deuteron_Phi");
      TBranch *Branch_Deuteron_TPC_Chi2overNDF	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_Chi2overNDF");
      TBranch *Branch_Deuteron_TPC_dEdx		    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_dEdx");
      TBranch *Branch_Deuteron_TPC_dEdx_Sigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_dEdx_Sigma");
      TBranch *Branch_Deuteron_TOF_Beta		    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Beta");
      TBranch *Branch_Deuteron_TOF_Beta_Sigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Beta_Sigma");
      TBranch *Branch_Deuteron_TOF_Mass2	    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Mass2");
      TBranch *Branch_Deuteron_TOF_Mass2_Sigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Mass2_Sigma");
      TBranch *Branch_Deuteron_ITS_dEdx		    = fTempTree_Deuteron->GetBranch("Deuteron_ITS_dEdx");
      TBranch *Branch_Deuteron_ITS_dEdx_Sigma	    = fTempTree_Deuteron->GetBranch("Deuteron_ITS_dEdx_Sigma");
      TBranch *Branch_Deuteron_DCAxy		    = fTempTree_Deuteron->GetBranch("Deuteron_DCAxy");
      TBranch *Branch_Deuteron_DCAz		    = fTempTree_Deuteron->GetBranch("Deuteron_DCAz");
      TBranch *Branch_Deuteron_Event_Centrality	    = fTempTree_Deuteron->GetBranch("Deuteron_Event_Centrality");
      TBranch *Branch_Deuteron_Event_PrimaryVertexZ = fTempTree_Deuteron->GetBranch("Deuteron_Event_PrimaryVertexZ");
      TBranch *Branch_Deuteron_TPC_nCrossedRows	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nCrossedRows");
      TBranch *Branch_Deuteron_TPC_nSharedCluster   = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nSharedCluster");
      TBranch *Branch_Deuteron_TPC_nClusterFindable = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nClusterFindable");
      TBranch *Branch_Deuteron_TPC_nCluster	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_nCluster");
      TBranch *Branch_Deuteron_ITS_nCluster	    = fTempTree_Deuteron->GetBranch("Deuteron_ITS_nCluster");
      TBranch *Branch_Deuteron_ID		    = fTempTree_Deuteron->GetBranch("Deuteron_ID");
      TBranch *Branch_Deuteron_Event_Identifier	    = fTempTree_Deuteron->GetBranch("Deuteron_Event_Identifier");


      Branch_Deuteron_px->SetAddress(&fDeuteron_px);
      Branch_Deuteron_py->SetAddress(&fDeuteron_py);
      Branch_Deuteron_pz->SetAddress(&fDeuteron_pz);
      Branch_Deuteron_pTPC->SetAddress(&fDeuteron_pTPC);
      Branch_Deuteron_Eta->SetAddress(&fDeuteron_Eta);
      Branch_Deuteron_Phi->SetAddress(&fDeuteron_Phi);
      Branch_Deuteron_TPC_Chi2overNDF->SetAddress(&fDeuteron_TPC_Chi2overNDF);
      Branch_Deuteron_TPC_dEdx->SetAddress(&fDeuteron_TPC_dEdx);
      Branch_Deuteron_TPC_dEdx_Sigma->SetAddress(&fDeuteron_TPC_dEdx_Sigma);
      Branch_Deuteron_TOF_Beta->SetAddress(&fDeuteron_TOF_Beta);
      Branch_Deuteron_TOF_Beta_Sigma->SetAddress(&fDeuteron_TOF_Beta_Sigma);
      Branch_Deuteron_TOF_Mass2->SetAddress(&fDeuteron_TOF_Mass2);
      Branch_Deuteron_TOF_Mass2_Sigma->SetAddress(&fDeuteron_TOF_Mass2_Sigma);
      Branch_Deuteron_ITS_dEdx->SetAddress(&fDeuteron_ITS_dEdx);
      Branch_Deuteron_ITS_dEdx_Sigma->SetAddress(&fDeuteron_ITS_dEdx_Sigma);
      Branch_Deuteron_DCAxy->SetAddress(&fDeuteron_DCAxy);
      Branch_Deuteron_DCAz->SetAddress(&fDeuteron_DCAz);
      Branch_Deuteron_Event_Centrality->SetAddress(&fDeuteron_Event_Centrality);
      Branch_Deuteron_Event_PrimaryVertexZ->SetAddress(&fDeuteron_Event_PrimaryVertexZ);
      Branch_Deuteron_TPC_nCrossedRows->SetAddress(&fDeuteron_TPC_nCrossedRows);
      Branch_Deuteron_TPC_nSharedCluster->SetAddress(&fDeuteron_TPC_nSharedCluster);
      Branch_Deuteron_TPC_nClusterFindable->SetAddress(&fDeuteron_TPC_nClusterFindable);
      Branch_Deuteron_TPC_nCluster->SetAddress(&fDeuteron_TPC_nCluster);
      Branch_Deuteron_ITS_nCluster->SetAddress(&fDeuteron_ITS_nCluster);
      Branch_Deuteron_ID->SetAddress(&fDeuteron_ID);
      Branch_Deuteron_Event_Identifier->SetAddress(&fDeuteron_Event_Identifier);


      Branch_Deuteron_px->SetAutoDelete(true);
      Branch_Deuteron_py->SetAutoDelete(true);
      Branch_Deuteron_pz->SetAutoDelete(true);
      Branch_Deuteron_pTPC->SetAutoDelete(true);
      Branch_Deuteron_Eta->SetAutoDelete(true);
      Branch_Deuteron_Phi->SetAutoDelete(true);
      Branch_Deuteron_TPC_Chi2overNDF->SetAutoDelete(true);
      Branch_Deuteron_TPC_dEdx->SetAutoDelete(true);
      Branch_Deuteron_TPC_dEdx_Sigma->SetAutoDelete(true);
      Branch_Deuteron_TOF_Beta->SetAutoDelete(true);
      Branch_Deuteron_TOF_Beta_Sigma->SetAutoDelete(true);
      Branch_Deuteron_TOF_Mass2->SetAutoDelete(true);
      Branch_Deuteron_TOF_Mass2_Sigma->SetAutoDelete(true);
      Branch_Deuteron_ITS_dEdx->SetAutoDelete(true);
      Branch_Deuteron_ITS_dEdx_Sigma->SetAutoDelete(true);
      Branch_Deuteron_DCAxy->SetAutoDelete(true);
      Branch_Deuteron_DCAz->SetAutoDelete(true);
      Branch_Deuteron_Event_Centrality->SetAutoDelete(true);
      Branch_Deuteron_Event_PrimaryVertexZ->SetAutoDelete(true);
      Branch_Deuteron_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_Deuteron_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_Deuteron_TPC_nClusterFindable->SetAutoDelete(true);
      Branch_Deuteron_TPC_nCluster->SetAutoDelete(true);
      Branch_Deuteron_ITS_nCluster->SetAutoDelete(true);
      Branch_Deuteron_ID->SetAutoDelete(true);
      Branch_Deuteron_Event_Identifier->SetAutoDelete(true);

      Branch_Deuteron_px->GetEntry(Deuteron);
      Branch_Deuteron_py->GetEntry(Deuteron);
      Branch_Deuteron_pz->GetEntry(Deuteron);
      Branch_Deuteron_pTPC->GetEntry(Deuteron);
      Branch_Deuteron_Eta->GetEntry(Deuteron);
      Branch_Deuteron_Phi->GetEntry(Deuteron);
      Branch_Deuteron_TPC_Chi2overNDF->GetEntry(Deuteron);
      Branch_Deuteron_TPC_dEdx->GetEntry(Deuteron);
      Branch_Deuteron_TPC_dEdx_Sigma->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Beta->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Beta_Sigma->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Mass2->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Mass2_Sigma->GetEntry(Deuteron);
      Branch_Deuteron_ITS_dEdx->GetEntry(Deuteron);
      Branch_Deuteron_ITS_dEdx_Sigma->GetEntry(Deuteron);
      Branch_Deuteron_DCAxy->GetEntry(Deuteron);
      Branch_Deuteron_DCAz->GetEntry(Deuteron);
      Branch_Deuteron_Event_Centrality->GetEntry(Deuteron);
      Branch_Deuteron_Event_PrimaryVertexZ->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nCrossedRows->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nSharedCluster->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nClusterFindable->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nCluster->GetEntry(Deuteron);
      Branch_Deuteron_ITS_nCluster->GetEntry(Deuteron);
      Branch_Deuteron_ID->GetEntry(Deuteron);
      Branch_Deuteron_Event_Identifier->GetEntry(Deuteron);
      fDeuteron_nParticlesPerEvent = nDeuteronsSelected;

      fSaveTree_Deuteron->Fill();


    } // end of loop (copy deuterons)


  } // end of same-event

  fTempTree_Proton->Delete();
  fTempTree_Deuteron->Delete();








  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antiproton selection loop +++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TTree     *fTempTree_AntiProton;
  float     AntiProton_px;
  float     AntiProton_py;
  float     AntiProton_pz;
  float     AntiProton_pTPC;
  float     AntiProton_Eta;
  float     AntiProton_Phi;
  float     AntiProton_TPC_Chi2overNDF;
  float     AntiProton_TPC_dEdx;
  float     AntiProton_TPC_dEdx_Sigma;
  float     AntiProton_TOF_Beta;
  float     AntiProton_TOF_Beta_Sigma;
  float     AntiProton_TOF_Mass2;
  float     AntiProton_TOF_Mass2_Sigma;
  float     AntiProton_ITS_dEdx;
  float     AntiProton_ITS_dEdx_Sigma;
  float     AntiProton_DCAxy;
  float     AntiProton_DCAz;
  float     AntiProton_Event_Centrality;
  float     AntiProton_Event_PrimaryVertexZ;
  unsigned short    AntiProton_TPC_nCrossedRows;
  unsigned short    AntiProton_TPC_nSharedCluster;
  unsigned short    AntiProton_TPC_nClusterFindable;
  unsigned short    AntiProton_TPC_nCluster;
  unsigned short    AntiProton_ITS_nCluster;
  unsigned int      AntiProton_ID;
  unsigned int      AntiProton_Event_Identifier;

  fTempTree_AntiProton = new TTree("fTempTree_AntiProton","fTempTree_AntiProton");
  fTempTree_AntiProton->Branch("AntiProton_px",&AntiProton_px,"AntiProton_px/f");
  fTempTree_AntiProton->Branch("AntiProton_py",&AntiProton_py,"AntiProton_py/f");
  fTempTree_AntiProton->Branch("AntiProton_pz",&AntiProton_pz,"AntiProton_pz/f");
  fTempTree_AntiProton->Branch("AntiProton_pTPC",&AntiProton_pTPC,"AntiProton_pTPC/f");
  fTempTree_AntiProton->Branch("AntiProton_Eta",&AntiProton_Eta,"AntiProton_Eta/f");
  fTempTree_AntiProton->Branch("AntiProton_Phi",&AntiProton_Phi,"AntiProton_Phi/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_Chi2overNDF",&AntiProton_TPC_Chi2overNDF,"AntiProton_TPC_Chi2overNDF/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_dEdx",&AntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_dEdx_Sigma",&AntiProton_TPC_dEdx_Sigma,"AntiProton_TPC_dEdx_Sigma/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Beta",&AntiProton_TOF_Beta,"AntiProton_TOF_Beta/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Beta_Sigma",&AntiProton_TOF_Beta_Sigma,"AntiProton_TOF_Beta_Sigma/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Mass2",&AntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Mass2_Sigma",&AntiProton_TOF_Mass2_Sigma,"AntiProton_TOF_Mass2_Sigma/f");
  fTempTree_AntiProton->Branch("AntiProton_ITS_dEdx",&AntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/f");
  fTempTree_AntiProton->Branch("AntiProton_ITS_dEdx_Sigma",&AntiProton_ITS_dEdx_Sigma,"AntiProton_ITS_dEdx_Sigma/f");
  fTempTree_AntiProton->Branch("AntiProton_DCAxy",&AntiProton_DCAxy,"AntiProton_DCAxy/f");
  fTempTree_AntiProton->Branch("AntiProton_DCAz",&AntiProton_DCAz,"AntiProton_DCAz/f");
  fTempTree_AntiProton->Branch("AntiProton_Event_Centrality",&AntiProton_Event_Centrality,"AntiProton_Event_Centrality/f");
  fTempTree_AntiProton->Branch("AntiProton_Event_PrimaryVertexZ",&AntiProton_Event_PrimaryVertexZ,"AntiProton_Event_PrimaryVertexZ/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&AntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&AntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nClusterFindable",&AntiProton_TPC_nClusterFindable,"AntiProton_TPC_nClusterFindable/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nCluster",&AntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_ITS_nCluster",&AntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_ID",&AntiProton_ID,"AntiProton_ID/i");
  fTempTree_AntiProton->Branch("AntiProton_Event_Identifier",&AntiProton_Event_Identifier,"AntiProton_Event_Identifier/i");

  unsigned short nAntiProtonsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// antiproton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antiproton cuts
    bool PassedAntiProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,false);
    if(!PassedAntiProtonCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    AntiProton_px		      = Track->Px();
    AntiProton_py		      = Track->Py();
    AntiProton_pz		      = Track->Pz();
    AntiProton_pTPC		      = Track->GetTPCmomentum();
    AntiProton_Eta		      = Track->Eta();
    AntiProton_Phi		      = Track->Phi();
    AntiProton_TPC_Chi2overNDF	      = Track->GetTPCchi2perNDF();
    AntiProton_TPC_dEdx		      = Track->GetTPCsignal();
    AntiProton_TPC_dEdx_Sigma	      = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    AntiProton_TOF_Beta		      = Track->GetTOFsignal();
    AntiProton_TOF_Beta_Sigma	      = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
    AntiProton_TOF_Mass2	      = CalculateMassSquareTOF(*Track);
    AntiProton_ITS_dEdx		      = Track->GetITSsignal();
    AntiProton_ITS_dEdx_Sigma	      = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kProton);
    AntiProton_DCAxy		      = DCAxy;
    AntiProton_DCAz		      = DCAz;
    AntiProton_Event_Centrality	      = Centrality;
    AntiProton_Event_PrimaryVertexZ   = PrimaryVertexZ;
    AntiProton_TPC_nCrossedRows	      = Track->GetTPCCrossedRows();
    AntiProton_TPC_nSharedCluster     = Track->GetTPCnclsS();
    AntiProton_TPC_nClusterFindable   = Track->GetTPCNclsF();
    AntiProton_TPC_nCluster	      = Track->GetTPCNcls();
    AntiProton_ITS_nCluster	      = Track->GetITSNcls();
    AntiProton_ID		      = Track->GetID();
    AntiProton_Event_Identifier	      = BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    fTempTree_AntiProton->Fill();
    nAntiProtonsSelected++;

  } // end of antiproton loop




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antideuteron selection loop +++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    TTree     *fTempTree_AntiDeuteron;
    float     AntiDeuteron_px;
    float     AntiDeuteron_py;
    float     AntiDeuteron_pz;
    float     AntiDeuteron_pTPC;
    float     AntiDeuteron_Eta;
    float     AntiDeuteron_Phi;
    float     AntiDeuteron_TPC_Chi2overNDF;
    float     AntiDeuteron_TPC_dEdx;
    float     AntiDeuteron_TPC_dEdx_Sigma;
    float     AntiDeuteron_TOF_Beta;
    float     AntiDeuteron_TOF_Beta_Sigma;
    float     AntiDeuteron_TOF_Mass2;
    float     AntiDeuteron_TOF_Mass2_Sigma;
    float     AntiDeuteron_ITS_dEdx;
    float     AntiDeuteron_ITS_dEdx_Sigma;
    float     AntiDeuteron_DCAxy;
    float     AntiDeuteron_DCAz;
    float     AntiDeuteron_Event_Centrality;
    float     AntiDeuteron_Event_PrimaryVertexZ;
    unsigned short  AntiDeuteron_TPC_nCrossedRows;
    unsigned short  AntiDeuteron_TPC_nSharedCluster;
    unsigned short  AntiDeuteron_TPC_nClusterFindable;
    unsigned short  AntiDeuteron_TPC_nCluster;
    unsigned short  AntiDeuteron_ITS_nCluster;
    unsigned int    AntiDeuteron_ID;
    unsigned int    AntiDeuteron_Event_Identifier;

  fTempTree_AntiDeuteron = new TTree("fTempTree_AntiDeuteron","fTempTree_AntiDeuteron");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_px",&AntiDeuteron_px,"AntiDeuteron_px/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_py",&AntiDeuteron_py,"AntiDeuteron_py/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pz",&AntiDeuteron_pz,"AntiDeuteron_pz/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&AntiDeuteron_pTPC,"AntiDeuteron_pTPC/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&AntiDeuteron_Eta,"AntiDeuteron_Eta/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&AntiDeuteron_Phi,"AntiDeuteron_Phi/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2overNDF",&AntiDeuteron_TPC_Chi2overNDF,"AntiDeuteron_TPC_Chi2overNDF/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&AntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_Sigma",&AntiDeuteron_TPC_dEdx_Sigma,"AntiDeuteron_TPC_dEdx_Sigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta",&AntiDeuteron_TOF_Beta,"AntiDeuteron_TOF_Beta/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta_Sigma",&AntiDeuteron_TOF_Beta_Sigma,"AntiDeuteron_TOF_Beta_Sigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&AntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_Sigma",&AntiDeuteron_TOF_Mass2_Sigma,"AntiDeuteron_TOF_Mass2_Sigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&AntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_Sigma",&AntiDeuteron_ITS_dEdx_Sigma,"AntiDeuteron_ITS_dEdx_Sigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&AntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&AntiDeuteron_DCAz,"AntiDeuteron_DCAz/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Event_Centrality",&AntiDeuteron_Event_Centrality,"AntiDeuteron_Event_Centrality/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Event_PrimaryVertexZ",&AntiDeuteron_Event_PrimaryVertexZ,"AntiDeuteron_Event_PrimaryVertexZ/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&AntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&AntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nClusterFindable",&AntiDeuteron_TPC_nClusterFindable,"AntiDeuteron_TPC_nClusterFindable/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&AntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&AntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ID",&AntiDeuteron_ID,"AntiDeuteron_ID/i");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&AntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/i");




  unsigned short nAntiDeuteronsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// antideuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antideuteron cuts
    bool PassedAntiDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,false);
    if(!PassedAntiDeuteronCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    AntiDeuteron_px			= Track->Px();
    AntiDeuteron_py			= Track->Py();
    AntiDeuteron_pz			= Track->Pz();
    AntiDeuteron_pTPC			= Track->GetTPCmomentum();
    AntiDeuteron_Eta			= Track->Eta();
    AntiDeuteron_Phi			= Track->Phi();
    AntiDeuteron_TPC_Chi2overNDF	= Track->GetTPCchi2perNDF();
    AntiDeuteron_TPC_dEdx		= Track->GetTPCsignal();
    AntiDeuteron_TPC_dEdx_Sigma		= fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    AntiDeuteron_TOF_Beta		= Track->GetTOFsignal();
    AntiDeuteron_TOF_Beta_Sigma		= fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
    AntiDeuteron_TOF_Mass2		= CalculateMassSquareTOF(*Track);
    AntiDeuteron_TOF_Mass2_Sigma	= CalculateDeuteronSigmaMassSquareTOF(Track->Pt(),AntiDeuteron_TOF_Mass2,false);
    AntiDeuteron_ITS_dEdx		= Track->GetITSsignal();
    AntiDeuteron_ITS_dEdx_Sigma		= fPIDResponse->NumberOfSigmasITS(Track,AliPID::kDeuteron);
    AntiDeuteron_DCAxy			= DCAxy;
    AntiDeuteron_DCAz			= DCAz;
    AntiDeuteron_Event_Centrality	= Centrality;
    AntiDeuteron_Event_PrimaryVertexZ	= PrimaryVertexZ;
    AntiDeuteron_TPC_nCrossedRows	= Track->GetTPCCrossedRows();
    AntiDeuteron_TPC_nSharedCluster	= Track->GetTPCnclsS();
    AntiDeuteron_TPC_nClusterFindable	= Track->GetTPCNclsF();
    AntiDeuteron_TPC_nCluster		= Track->GetTPCNcls();
    AntiDeuteron_ITS_nCluster		= Track->GetITSNcls();
    AntiDeuteron_ID			= Track->GetID();
    AntiDeuteron_Event_Identifier	= BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    fTempTree_AntiDeuteron->Fill();
    nAntiDeuteronsSelected++;

  } // end of antideuteron loop




  if((nAntiProtonsSelected > 0) && (nAntiDeuteronsSelected > 0)){

    for(int AntiProton = 0; AntiProton < nAntiProtonsSelected; AntiProton++){

      TBranch *Branch_AntiProton_px			  = fTempTree_AntiProton->GetBranch("AntiProton_px");
      TBranch *Branch_AntiProton_py			  = fTempTree_AntiProton->GetBranch("AntiProton_py");
      TBranch *Branch_AntiProton_pz			  = fTempTree_AntiProton->GetBranch("AntiProton_pz");
      TBranch *Branch_AntiProton_pTPC		  = fTempTree_AntiProton->GetBranch("AntiProton_pTPC");
      TBranch *Branch_AntiProton_Eta		  = fTempTree_AntiProton->GetBranch("AntiProton_Eta");
      TBranch *Branch_AntiProton_Phi		  = fTempTree_AntiProton->GetBranch("AntiProton_Phi");
      TBranch *Branch_AntiProton_TPC_Chi2overNDF	  = fTempTree_AntiProton->GetBranch("AntiProton_TPC_Chi2overNDF");
      TBranch *Branch_AntiProton_TPC_dEdx		  = fTempTree_AntiProton->GetBranch("AntiProton_TPC_dEdx");
      TBranch *Branch_AntiProton_TPC_dEdx_Sigma	  = fTempTree_AntiProton->GetBranch("AntiProton_TPC_dEdx_Sigma");
      TBranch *Branch_AntiProton_TOF_Beta		  = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Beta");
      TBranch *Branch_AntiProton_TOF_Beta_Sigma	  = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Beta_Sigma");
      TBranch *Branch_AntiProton_TOF_Mass2		  = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Mass2");
      TBranch *Branch_AntiProton_TOF_Mass2_Sigma	  = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Mass2_Sigma");
      TBranch *Branch_AntiProton_ITS_dEdx		  = fTempTree_AntiProton->GetBranch("AntiProton_ITS_dEdx");
      TBranch *Branch_AntiProton_ITS_dEdx_Sigma	  = fTempTree_AntiProton->GetBranch("AntiProton_ITS_dEdx_Sigma");
      TBranch *Branch_AntiProton_DCAxy		  = fTempTree_AntiProton->GetBranch("AntiProton_DCAxy");
      TBranch *Branch_AntiProton_DCAz		  = fTempTree_AntiProton->GetBranch("AntiProton_DCAz");
      TBranch *Branch_AntiProton_Event_Centrality	  = fTempTree_AntiProton->GetBranch("AntiProton_Event_Centrality");
      TBranch *Branch_AntiProton_Event_PrimaryVertexZ = fTempTree_AntiProton->GetBranch("AntiProton_Event_PrimaryVertexZ");
      TBranch *Branch_AntiProton_TPC_nCrossedRows	  = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nCrossedRows");
      TBranch *Branch_AntiProton_TPC_nSharedCluster	  = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nSharedCluster");
      TBranch *Branch_AntiProton_TPC_nClusterFindable = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nClusterFindable");
      TBranch *Branch_AntiProton_TPC_nCluster	  = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nCluster");
      TBranch *Branch_AntiProton_ITS_nCluster	  = fTempTree_AntiProton->GetBranch("AntiProton_ITS_nCluster");
      TBranch *Branch_AntiProton_ID			  = fTempTree_AntiProton->GetBranch("AntiProton_ID");
      TBranch *Branch_AntiProton_Event_Identifier	  = fTempTree_AntiProton->GetBranch("AntiProton_Event_Identifier");


      Branch_AntiProton_px->SetAddress(&fAntiProton_px);
      Branch_AntiProton_py->SetAddress(&fAntiProton_py);
      Branch_AntiProton_pz->SetAddress(&fAntiProton_pz);
      Branch_AntiProton_pTPC->SetAddress(&fAntiProton_pTPC);
      Branch_AntiProton_Eta->SetAddress(&fAntiProton_Eta);
      Branch_AntiProton_Phi->SetAddress(&fAntiProton_Phi);
      Branch_AntiProton_TPC_Chi2overNDF->SetAddress(&fAntiProton_TPC_Chi2overNDF);
      Branch_AntiProton_TPC_dEdx->SetAddress(&fAntiProton_TPC_dEdx);
      Branch_AntiProton_TPC_dEdx_Sigma->SetAddress(&fAntiProton_TPC_dEdx_Sigma);
      Branch_AntiProton_TOF_Beta->SetAddress(&fAntiProton_TOF_Beta);
      Branch_AntiProton_TOF_Beta_Sigma->SetAddress(&fAntiProton_TOF_Beta_Sigma);
      Branch_AntiProton_TOF_Mass2->SetAddress(&fAntiProton_TOF_Mass2);
      Branch_AntiProton_TOF_Mass2_Sigma->SetAddress(&fAntiProton_TOF_Mass2_Sigma);
      Branch_AntiProton_ITS_dEdx->SetAddress(&fAntiProton_ITS_dEdx);
      Branch_AntiProton_ITS_dEdx_Sigma->SetAddress(&fAntiProton_ITS_dEdx_Sigma);
      Branch_AntiProton_DCAxy->SetAddress(&fAntiProton_DCAxy);
      Branch_AntiProton_DCAz->SetAddress(&fAntiProton_DCAz);
      Branch_AntiProton_Event_Centrality->SetAddress(&fAntiProton_Event_Centrality);
      Branch_AntiProton_Event_PrimaryVertexZ->SetAddress(&fAntiProton_Event_PrimaryVertexZ);
      Branch_AntiProton_TPC_nCrossedRows->SetAddress(&fAntiProton_TPC_nCrossedRows);
      Branch_AntiProton_TPC_nSharedCluster->SetAddress(&fAntiProton_TPC_nSharedCluster);
      Branch_AntiProton_TPC_nClusterFindable->SetAddress(&fAntiProton_TPC_nClusterFindable);
      Branch_AntiProton_TPC_nCluster->SetAddress(&fAntiProton_TPC_nCluster);
      Branch_AntiProton_ITS_nCluster->SetAddress(&fAntiProton_ITS_nCluster);
      Branch_AntiProton_ID->SetAddress(&fAntiProton_ID);
      Branch_AntiProton_Event_Identifier->SetAddress(&fAntiProton_Event_Identifier);


      Branch_AntiProton_px->SetAutoDelete(true);
      Branch_AntiProton_py->SetAutoDelete(true);
      Branch_AntiProton_pz->SetAutoDelete(true);
      Branch_AntiProton_pTPC->SetAutoDelete(true);
      Branch_AntiProton_Eta->SetAutoDelete(true);
      Branch_AntiProton_Phi->SetAutoDelete(true);
      Branch_AntiProton_TPC_Chi2overNDF->SetAutoDelete(true);
      Branch_AntiProton_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiProton_TPC_dEdx_Sigma->SetAutoDelete(true);
      Branch_AntiProton_TOF_Beta->SetAutoDelete(true);
      Branch_AntiProton_TOF_Beta_Sigma->SetAutoDelete(true);
      Branch_AntiProton_TOF_Mass2->SetAutoDelete(true);
      Branch_AntiProton_TOF_Mass2_Sigma->SetAutoDelete(true);
      Branch_AntiProton_ITS_dEdx->SetAutoDelete(true);
      Branch_AntiProton_ITS_dEdx_Sigma->SetAutoDelete(true);
      Branch_AntiProton_DCAxy->SetAutoDelete(true);
      Branch_AntiProton_DCAz->SetAutoDelete(true);
      Branch_AntiProton_Event_Centrality->SetAutoDelete(true);
      Branch_AntiProton_Event_PrimaryVertexZ->SetAutoDelete(true);
      Branch_AntiProton_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_AntiProton_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_AntiProton_TPC_nClusterFindable->SetAutoDelete(true);
      Branch_AntiProton_TPC_nCluster->SetAutoDelete(true);
      Branch_AntiProton_ITS_nCluster->SetAutoDelete(true);
      Branch_AntiProton_ID->SetAutoDelete(true);
      Branch_AntiProton_Event_Identifier->SetAutoDelete(true);

      Branch_AntiProton_px->GetEntry(AntiProton);
      Branch_AntiProton_py->GetEntry(AntiProton);
      Branch_AntiProton_pz->GetEntry(AntiProton);
      Branch_AntiProton_pTPC->GetEntry(AntiProton);
      Branch_AntiProton_Eta->GetEntry(AntiProton);
      Branch_AntiProton_Phi->GetEntry(AntiProton);
      Branch_AntiProton_TPC_Chi2overNDF->GetEntry(AntiProton);
      Branch_AntiProton_TPC_dEdx->GetEntry(AntiProton);
      Branch_AntiProton_TPC_dEdx_Sigma->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Beta->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Beta_Sigma->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Mass2->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Mass2_Sigma->GetEntry(AntiProton);
      Branch_AntiProton_ITS_dEdx->GetEntry(AntiProton);
      Branch_AntiProton_ITS_dEdx_Sigma->GetEntry(AntiProton);
      Branch_AntiProton_DCAxy->GetEntry(AntiProton);
      Branch_AntiProton_DCAz->GetEntry(AntiProton);
      Branch_AntiProton_Event_Centrality->GetEntry(AntiProton);
      Branch_AntiProton_Event_PrimaryVertexZ->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nCrossedRows->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nSharedCluster->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nClusterFindable->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nCluster->GetEntry(AntiProton);
      Branch_AntiProton_ITS_nCluster->GetEntry(AntiProton);
      Branch_AntiProton_ID->GetEntry(AntiProton);
      Branch_AntiProton_Event_Identifier->GetEntry(AntiProton);
      fAntiProton_nParticlesPerEvent = nAntiProtonsSelected;

      fSaveTree_AntiProton->Fill();

    } // end of loop (copy antiprotons)


    for(int AntiDeuteron = 0; AntiDeuteron < nAntiDeuteronsSelected; AntiDeuteron++){

      TBranch *Branch_AntiDeuteron_px			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_px");
      TBranch *Branch_AntiDeuteron_py			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_py");
      TBranch *Branch_AntiDeuteron_pz			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pz");
      TBranch *Branch_AntiDeuteron_pTPC			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pTPC");
      TBranch *Branch_AntiDeuteron_Eta			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Eta");
      TBranch *Branch_AntiDeuteron_Phi			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Phi");
      TBranch *Branch_AntiDeuteron_TPC_Chi2overNDF	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_Chi2overNDF");
      TBranch *Branch_AntiDeuteron_TPC_dEdx		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_dEdx");
      TBranch *Branch_AntiDeuteron_TPC_dEdx_Sigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_dEdx_Sigma");
      TBranch *Branch_AntiDeuteron_TOF_Beta		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Beta");
      TBranch *Branch_AntiDeuteron_TOF_Beta_Sigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Beta_Sigma");
      TBranch *Branch_AntiDeuteron_TOF_Mass2		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Mass2");
      TBranch *Branch_AntiDeuteron_TOF_Mass2_Sigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Mass2_Sigma");
      TBranch *Branch_AntiDeuteron_ITS_dEdx		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ITS_dEdx");
      TBranch *Branch_AntiDeuteron_ITS_dEdx_Sigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ITS_dEdx_Sigma");
      TBranch *Branch_AntiDeuteron_DCAxy		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_DCAxy");
      TBranch *Branch_AntiDeuteron_DCAz			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_DCAz");
      TBranch *Branch_AntiDeuteron_Event_Centrality	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Event_Centrality");
      TBranch *Branch_AntiDeuteron_Event_PrimaryVertexZ	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Event_PrimaryVertexZ");
      TBranch *Branch_AntiDeuteron_TPC_nCrossedRows	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nCrossedRows");
      TBranch *Branch_AntiDeuteron_TPC_nSharedCluster   = fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nSharedCluster");
      TBranch *Branch_AntiDeuteron_TPC_nClusterFindable = fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nClusterFindable");
      TBranch *Branch_AntiDeuteron_TPC_nCluster		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_nCluster");
      TBranch *Branch_AntiDeuteron_ITS_nCluster		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ITS_nCluster");
      TBranch *Branch_AntiDeuteron_ID			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ID");
      TBranch *Branch_AntiDeuteron_Event_Identifier	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Event_Identifier");


      Branch_AntiDeuteron_px->SetAddress(&fAntiDeuteron_px);
      Branch_AntiDeuteron_py->SetAddress(&fAntiDeuteron_py);
      Branch_AntiDeuteron_pz->SetAddress(&fAntiDeuteron_pz);
      Branch_AntiDeuteron_pTPC->SetAddress(&fAntiDeuteron_pTPC);
      Branch_AntiDeuteron_Eta->SetAddress(&fAntiDeuteron_Eta);
      Branch_AntiDeuteron_Phi->SetAddress(&fAntiDeuteron_Phi);
      Branch_AntiDeuteron_TPC_Chi2overNDF->SetAddress(&fAntiDeuteron_TPC_Chi2overNDF);
      Branch_AntiDeuteron_TPC_dEdx->SetAddress(&fAntiDeuteron_TPC_dEdx);
      Branch_AntiDeuteron_TPC_dEdx_Sigma->SetAddress(&fAntiDeuteron_TPC_dEdx_Sigma);
      Branch_AntiDeuteron_TOF_Beta->SetAddress(&fAntiDeuteron_TOF_Beta);
      Branch_AntiDeuteron_TOF_Beta_Sigma->SetAddress(&fAntiDeuteron_TOF_Beta_Sigma);
      Branch_AntiDeuteron_TOF_Mass2->SetAddress(&fAntiDeuteron_TOF_Mass2);
      Branch_AntiDeuteron_TOF_Mass2_Sigma->SetAddress(&fAntiDeuteron_TOF_Mass2_Sigma);
      Branch_AntiDeuteron_ITS_dEdx->SetAddress(&fAntiDeuteron_ITS_dEdx);
      Branch_AntiDeuteron_ITS_dEdx_Sigma->SetAddress(&fAntiDeuteron_ITS_dEdx_Sigma);
      Branch_AntiDeuteron_DCAxy->SetAddress(&fAntiDeuteron_DCAxy);
      Branch_AntiDeuteron_DCAz->SetAddress(&fAntiDeuteron_DCAz);
      Branch_AntiDeuteron_Event_Centrality->SetAddress(&fAntiDeuteron_Event_Centrality);
      Branch_AntiDeuteron_Event_PrimaryVertexZ->SetAddress(&fAntiDeuteron_Event_PrimaryVertexZ);
      Branch_AntiDeuteron_TPC_nCrossedRows->SetAddress(&fAntiDeuteron_TPC_nCrossedRows);
      Branch_AntiDeuteron_TPC_nSharedCluster->SetAddress(&fAntiDeuteron_TPC_nSharedCluster);
      Branch_AntiDeuteron_TPC_nClusterFindable->SetAddress(&fAntiDeuteron_TPC_nClusterFindable);
      Branch_AntiDeuteron_TPC_nCluster->SetAddress(&fAntiDeuteron_TPC_nCluster);
      Branch_AntiDeuteron_ITS_nCluster->SetAddress(&fAntiDeuteron_ITS_nCluster);
      Branch_AntiDeuteron_ID->SetAddress(&fAntiDeuteron_ID);
      Branch_AntiDeuteron_Event_Identifier->SetAddress(&fAntiDeuteron_Event_Identifier);


      Branch_AntiDeuteron_px->SetAutoDelete(true);
      Branch_AntiDeuteron_py->SetAutoDelete(true);
      Branch_AntiDeuteron_pz->SetAutoDelete(true);
      Branch_AntiDeuteron_pTPC->SetAutoDelete(true);
      Branch_AntiDeuteron_Eta->SetAutoDelete(true);
      Branch_AntiDeuteron_Phi->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_Chi2overNDF->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_dEdx_Sigma->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Beta->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Beta_Sigma->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Mass2->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Mass2_Sigma->SetAutoDelete(true);
      Branch_AntiDeuteron_ITS_dEdx->SetAutoDelete(true);
      Branch_AntiDeuteron_ITS_dEdx_Sigma->SetAutoDelete(true);
      Branch_AntiDeuteron_DCAxy->SetAutoDelete(true);
      Branch_AntiDeuteron_DCAz->SetAutoDelete(true);
      Branch_AntiDeuteron_Event_Centrality->SetAutoDelete(true);
      Branch_AntiDeuteron_Event_PrimaryVertexZ->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nClusterFindable->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_nCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_ITS_nCluster->SetAutoDelete(true);
      Branch_AntiDeuteron_ID->SetAutoDelete(true);
      Branch_AntiDeuteron_Event_Identifier->SetAutoDelete(true);

      Branch_AntiDeuteron_px->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_py->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_pz->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_pTPC->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Eta->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Phi->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_Chi2overNDF->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_dEdx->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_dEdx_Sigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Beta->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Beta_Sigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Mass2->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Mass2_Sigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ITS_dEdx->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ITS_dEdx_Sigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_DCAxy->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_DCAz->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Event_Centrality->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Event_PrimaryVertexZ->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nCrossedRows->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nSharedCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nClusterFindable->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ITS_nCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ID->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Event_Identifier->GetEntry(AntiDeuteron);
      fAntiDeuteron_nParticlesPerEvent = nAntiDeuteronsSelected;

      fSaveTree_AntiDeuteron->Fill();


    } // end of loop (copy antideuterons)


  } // end of same-event

  fTempTree_AntiProton->Delete();
  fTempTree_AntiDeuteron->Delete();












  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiProton);
  PostData(4,fSaveTree_AntiDeuteron);

} // end of UserExec















void AliAnalysisTask_pd_CreateTrees_PairsOnly::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateBetaTOF(AliAODTrack &track)
{

  double length = track.GetIntegratedLength(); // cm
  if(length <= 350.0) return -999.0;

  double c = TMath::C(); // m/s
  double end_time = track.GetTOFsignal(); // ps
  double start_time = fPIDResponse->GetTOFResponse().GetStartTime(track.GetTPCmomentum()); // ps
  double time = (end_time - start_time) * 1e-12; // ps -> s

  double velocity = (length*0.01) / time; // m/s
  double beta = velocity / c;

  return beta;

} // end of CalculateBetaTOF






// calculate the mass² TOF
double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateMassSquareTOF(AliAODTrack &track)
{

  double p = track.P();
  double beta = CalculateBetaTOF(track);
  double mass2 = -999.0;

  if(beta > 0.0){

    mass2 = (1/(beta*beta)-1) * (p*p);

  }

  return mass2;

} // end of CalculateMassSquareTOF








double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateDeuteronSigmaMassSquareTOF(double pT, double massSq, bool isMatter)
{

  // for Pb-Pb (LHC18q and LHC18r)

  double ParameterMatter[35][2] = {
{4.17615,0.298251},
{4.08097,0.297396},
{3.95197,0.240426},
{3.82068,0.149232},
{3.74164,0.118375},
{3.69264,0.11247},
{3.65695,0.112983},
{3.63148,0.111041},
{3.61447,0.111148},
{3.60399,0.112195},
{3.59677,0.113264},
{3.59217,0.114995},
{3.59126,0.118023},
{3.59045,0.120715},
{3.58861,0.122859},
{3.58585,0.126467},
{3.58287,0.12996},
{3.5818,0.134162},
{3.58132,0.138716},
{3.58191,0.142358},
{3.58124,0.147029},
{3.5818,0.152781},
{3.58151,0.156665},
{3.58148,0.160992},
{3.58288,0.16659},
{3.58375,0.171199},
{3.58526,0.1769},
{3.58636,0.182282},
{3.58728,0.188404},
{3.58966,0.194788},
{3.59245,0.199839},
{3.59187,0.206191},
{3.59449,0.213257},
{3.59563,0.216757},
}; // end of ParameterMatter definition


  double ParameterAntiMatter[35][2] = {
{4.34701,0.293789},
{4.1721,0.270784},
{3.97304,0.210971},
{3.83358,0.16},
{3.74727,0.130793},
{3.69333,0.119528},
{3.65658,0.115693},
{3.63029,0.113118},
{3.6127,0.112189},
{3.60345,0.113028},
{3.59819,0.115984},
{3.59424,0.116708},
{3.59396,0.120272},
{3.59399,0.123617},
{3.59312,0.12665},
{3.59111,0.129108},
{3.58911,0.13153},
{3.58753,0.135724},
{3.58664,0.140496},
{3.58885,0.146367},
{3.58733,0.150429},
{3.58957,0.152747},
{3.58825,0.156932},
{3.59141,0.162311},
{3.59057,0.166413},
{3.59221,0.172248},
{3.59319,0.178037},
{3.59445,0.184472},
{3.5954,0.190762},
{3.59904,0.198087},
{3.6033,0.20339},
{3.60191,0.207736},
{3.60263,0.213331},
{3.60588,0.214263},
}; // end of ParameterAntiMatter definition

  int row = 0; // pt < 0.5         

  if(pT > 0.6) row = 1;
  if(pT > 0.7) row = 2;
  if(pT > 0.8) row = 3;
  if(pT > 0.9) row = 4;
  if(pT > 1.0) row = 5;
  if(pT > 1.1) row = 6;
  if(pT > 1.2) row = 7;
  if(pT > 1.3) row = 8;
  if(pT > 1.4) row = 9;
  if(pT > 1.5) row = 10;
  if(pT > 1.6) row = 11;
  if(pT > 1.7) row = 12;
  if(pT > 1.8) row = 13;
  if(pT > 1.9) row = 14;
  if(pT > 2.0) row = 15;
  if(pT > 2.1) row = 16;
  if(pT > 2.2) row = 17;
  if(pT > 2.3) row = 18;
  if(pT > 2.4) row = 19;
  if(pT > 2.5) row = 20;
  if(pT > 2.6) row = 21;
  if(pT > 2.7) row = 22;
  if(pT > 2.8) row = 23;
  if(pT > 2.9) row = 24;
  if(pT > 3.0) row = 25;
  if(pT > 3.1) row = 26;
  if(pT > 3.2) row = 27;
  if(pT > 3.3) row = 28;
  if(pT > 3.4) row = 29;
  if(pT > 3.5) row = 30;
  if(pT > 3.6) row = 31;
  if(pT > 3.7) row = 32;
  if(pT > 3.8) row = 33;
  if(pT > 3.9) row = 34;

  double mean = 0.0;
  double sigma = 0.0;

  if(isMatter) mean   = ParameterMatter[row][0];
  if(isMatter) sigma  = ParameterMatter[row][1];

  if(!isMatter) mean   = ParameterAntiMatter[row][0];
  if(!isMatter) sigma  = ParameterAntiMatter[row][1];

  if(TMath::Abs(sigma) < 0.00001) return -999.0;

  double nSigma = (massSq - mean)/(sigma);
  return nSigma;

} // end of CalculateDeuteronSigmaMassSquareTOF






// apply track cuts for protons and antiprotons
bool AliAnalysisTask_pd_CreateTrees_PairsOnly::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter)
{

  bool PassedParticleCuts = false;

  // define proton and antiproton track cuts
  double Proton_pT_min = 0.4; // GeV/c
  double Proton_pT_max = 4.1; // GeV/c

  double Proton_eta_min = -0.9;
  double Proton_eta_max = +0.9;

  double Proton_DCAxy_max = 0.2; // cm
  double Proton_DCAz_max = 0.3; // cm

  double Proton_TPC_RatioRowsCluster_min = 0.75;

  double Proton_TPC_nSigma_max = 4.0;
  double Proton_TOF_nSigma_max = 4.0;

  double Proton_TOF_nSigma_max_low_p = 7.0;

  int Proton_TPC_nCluster_min = 70;

  int Proton_TPC_nCrossedRows_min = 60;

  int Proton_TPC_nSharedCluster_max = 3;

  int Proton_FilterBit = 7; // FilterBit 7 = 128

  int Proton_ITS_nCluster_min = 1;
  double Proton_ITS_nSigma_max = 60;

  bool UseITS = true;
  double Proton_TPC_Threshold = 0.70;


  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;

  // apple TPC sigma cut
  double nSigmaTPC = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
  if(TMath::Abs(nSigmaTPC) > Proton_TPC_nSigma_max) return PassedParticleCuts;

  // get DCA infotmation
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  double DCAxy = xv[0];
  double DCAz = xv[1];

  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Proton_DCAxy_max) return PassedParticleCuts;

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Proton_DCAz_max) return PassedParticleCuts;

  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  double p    = Track.P();
  double pT   = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  // check if TOF information is available above threshold
  if((pTPC >= Proton_TPC_Threshold) && (!TOFisOK)) return PassedParticleCuts;

  double nSigmaTOF = -999.0;
  if(TOFisOK) nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kProton);

  // apply TOF low p cut
  if(TOFisOK)
  {

    if((pTPC < Proton_TPC_Threshold) && (TMath::Abs(nSigmaTOF) > Proton_TOF_nSigma_max_low_p)) return PassedParticleCuts;

  }

  // apply TOF high p cut
  if(TOFisOK)
  {

    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kProton);

    if((pTPC >= Proton_TPC_Threshold) && (TMath::Abs(nSigmaTOF) > Proton_TOF_nSigma_max)) return PassedParticleCuts;

  }

/*
  // reject tracks with better sigma for other particles
  if((pTPC >= Proton_TPC_Threshold) && (TOFisOK)){

    // apply TPC and TOF sigma cut
    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kProton);

    double nSigmaTPC_Electron = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kElectron);
    double nSigmaTPC_Deuteron = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
    double nSigmaTPC_Muon     = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kMuon);
    double nSigmaTPC_Kaon     = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kKaon);
    double nSigmaTPC_Pion     = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kPion);

    double nSigmaTOF_Electron = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kElectron);
    double nSigmaTOF_Deuteron = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);
    double nSigmaTOF_Muon     = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kMuon);
    double nSigmaTOF_Kaon     = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kKaon);
    double nSigmaTOF_Pion     = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kPion);

    // reject candidates which have a smaller sigma for other particles
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Electron)) && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Electron))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Deuteron)) && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Deuteron))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Muon))	  && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Muon))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Kaon))	  && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Kaon))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Pion))	  && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Pion))) return PassedParticleCuts;

  }
*/

  // apply FilterBit cut
  if(!Track.TestFilterBit(Proton_FilterBit)) return PassedParticleCuts;

  // apply pT cut
  if(pT < Proton_pT_min || pT > Proton_pT_max) return PassedParticleCuts;

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(eta < Proton_eta_min || eta > Proton_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  int nClusterTPC = Track.GetNcls(1);
  if(nClusterTPC < Proton_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  int nCrossedRowsTPC = Track.GetTPCCrossedRows();
  if(nCrossedRowsTPC < Proton_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply zero shared cluster cut for TPC
  int nSharedClusterTPC = Track.GetTPCnclsS();
  if(nSharedClusterTPC > Proton_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  int nClusterTPCfindable = Track.GetTPCNclsF();
  double RatioRowsFindableClusterTPC = -999.0;
  if(nClusterTPCfindable > 0) RatioRowsFindableClusterTPC = ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);
  if(RatioRowsFindableClusterTPC < Proton_TPC_RatioRowsCluster_min) return PassedParticleCuts;


  // apply ITS cluster cut
  double nClusterITS = Track.GetITSNcls();
  if((UseITS) && (nClusterITS <= Proton_ITS_nCluster_min)) return PassedParticleCuts;

  // apply ITS dEdx cut below proton band
  if((UseITS))
  {

    TF1 *fdEdxProtonITS = new TF1("fdEdxProtonITS","0.5*[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.);
    fdEdxProtonITS->SetParameters(2.36861e-07,-55831.1,-238672,9.55834,17081,1,0.93827208816);

    if(Track.GetITSsignal() < fdEdxProtonITS->Eval(p)) return PassedParticleCuts;

  }

  // apply ITS dEdx cut above proton band
  if((UseITS))
  {

    TF1 *fdEdxProtonITS = new TF1("fdEdxProtonITS","1.5*[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.);
    fdEdxProtonITS->SetParameters(2.36861e-07,-55831.1,-238672,9.55834,17081,1,0.93827208816);

    if(Track.GetITSsignal() > fdEdxProtonITS->Eval(p)) return PassedParticleCuts;

  }


  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckProtonCuts


// apply track cuts for deuterons and antideuterons
bool AliAnalysisTask_pd_CreateTrees_PairsOnly::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter)
{

  bool PassedParticleCuts = false;

  // define deuteron and antideuteron track cuts
  double Deuteron_pT_min = 0.4;
  double Deuteron_pT_max = 1.6;
  double Deuteron_eta_min = -0.9;
  double Deuteron_eta_max = +0.9;
  double Deuteron_DCAxy_max = 0.2; // cm
  double Deuteron_DCAz_max = 0.3; // cm

  double Deuteron_TPC_RatioRowsCluster_min = 0.7;
  double Deuteron_TPC_nSigma_max = 4.0;
  double Deuteron_TOF_nSigma_max = 4.0;

  double Deuteron_TOF_nSigma_max_low_p = 7.0;

  int Deuteron_TPC_nCluster_min = 70;
  int Deuteron_TPC_nCrossedRows_min = 60;
  int Deuteron_TPC_nSharedCluster_max = 3;
  int Deuteron_FilterBit = 7; // FilterBit 7 = 128
  int Deuteron_ITS_nCluster_min = 1;


  double Deuteron_TPC_Threshold = 1.0;

  double Deuteron_TOF_MassSquare_nSigma_max = 4.0;

  bool UseBetaTOF = false;
  bool UseMassSquareTOF = true;
  bool UseITS = true;


  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;

  // apply TPC nSigma cut
  double nSigmaTPC = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
  if(TMath::Abs(nSigmaTPC) > Deuteron_TPC_nSigma_max) return PassedParticleCuts;

  // get DCA information
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  double DCAxy = xv[0];
  double DCAz = xv[1];
  
  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Deuteron_DCAxy_max) return PassedParticleCuts;

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Deuteron_DCAz_max) return PassedParticleCuts;

  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  double p = Track.P();
  double pT = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  // check TOF status above threshold
  if((pTPC >= Deuteron_TPC_Threshold) && (!TOFisOK)) return PassedParticleCuts;

  double nSigmaTOF = -999.0;
    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);

  // apply TOF low p cut
  if(TOFisOK)
  {

    if((pTPC < Deuteron_TPC_Threshold) && (TMath::Abs(nSigmaTOF) > Deuteron_TOF_nSigma_max_low_p)) return PassedParticleCuts;

  }


  if((UseBetaTOF) && (pTPC >= Deuteron_TPC_Threshold))
  {

    // apply TPC and TOF sigma cut

    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);
    if((TMath::Abs(nSigmaTPC) > Deuteron_TPC_nSigma_max) || (TMath::Abs(nSigmaTOF) > Deuteron_TOF_nSigma_max)) return PassedParticleCuts;



    // reject tracks with better sigma for other particles
    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);

    double nSigmaTPC_Electron = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kElectron);
    double nSigmaTPC_Deuteron = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
    double nSigmaTPC_Muon     = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kMuon);
    double nSigmaTPC_Kaon     = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kKaon);
    double nSigmaTPC_Pion     = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kPion);

    double nSigmaTOF_Electron = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kElectron);
    double nSigmaTOF_Deuteron = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);
    double nSigmaTOF_Muon     = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kMuon);
    double nSigmaTOF_Kaon     = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kKaon);
    double nSigmaTOF_Pion     = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kPion);

    // reject candidates which have a smaller sigma for other particles
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Electron)) && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Electron))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Deuteron)) && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Deuteron))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Muon))	  && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Muon))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Kaon))	  && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Kaon))) return PassedParticleCuts;
    if((TMath::Abs(nSigmaTPC) > TMath::Abs(nSigmaTPC_Pion))	  && (TMath::Abs(nSigmaTOF) > TMath::Abs(nSigmaTOF_Pion))) return PassedParticleCuts;


  } // end of UseBetaTOF


  if(UseMassSquareTOF)
  {
  
    if(pTPC >= Deuteron_TPC_Threshold)
    {

    double massSq = CalculateMassSquareTOF(Track);
    double nSigmaTOFmSq = CalculateDeuteronSigmaMassSquareTOF(pT,massSq,isMatter);

    if(TMath::Abs(nSigmaTOFmSq) > Deuteron_TOF_MassSquare_nSigma_max) return PassedParticleCuts;

    } // end of UseMassSquareTOF


  }

  // apply FilterBit cut
  if(!Track.TestFilterBit(Deuteron_FilterBit)) return PassedParticleCuts;

  // apply pT cut
  if(pT < Deuteron_pT_min || pT > Deuteron_pT_max) return PassedParticleCuts;

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(eta < Deuteron_eta_min || eta > Deuteron_eta_max) return PassedParticleCuts;

  // apply cluster cut for TPC
  int nClusterTPC = Track.GetNcls(1);
  if(nClusterTPC < Deuteron_TPC_nCluster_min) return PassedParticleCuts;

  // apply crossed rows cut for TPC
  int nCrossedRowsTPC = Track.GetTPCCrossedRows();
  if(nCrossedRowsTPC < Deuteron_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply zero shared cluster cut for TPC
  int nSharedClusterTPC = Track.GetTPCnclsS();
  if(nSharedClusterTPC > Deuteron_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  int nClusterTPCfindable = Track.GetTPCNclsF();
  double RatioRowsFindableClusterTPC = -999.0;
  if(nClusterTPCfindable > 0) RatioRowsFindableClusterTPC = ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);
  if(RatioRowsFindableClusterTPC < Deuteron_TPC_RatioRowsCluster_min) return PassedParticleCuts;


  // check if ITS information is available
  //AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  //bool ITSisOK = false;
  //if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;

  // apply ITS cluster cut
  if(UseITS)
  {

    double nClusterITS = Track.GetITSNcls();
    if(nClusterITS <= Deuteron_ITS_nCluster_min) return PassedParticleCuts;

  }


  // apply ITS dEdx cut
  if(UseITS)
  {

    TF1 *fdEdxDeuteronITS = new TF1("fdEdxDeuteronITS","0.7*[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.);
    fdEdxDeuteronITS->SetParameters(7.41722e-06,-55831.1,-238672,11249.3,19828.9,1,1.8756129425);

    if(Track.GetITSsignal() < fdEdxDeuteronITS->Eval(p)) return PassedParticleCuts;

  }



  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckDeuteronCuts



