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

#include "AliAnalysisTask_pd.h"

using namespace std;
ClassImp(AliAnalysisTask_pd) 






AliAnalysisTask_pd::AliAnalysisTask_pd() : AliAnalysisTaskSE(),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(0),
  fTree_Proton(0),
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
  fProton_Event_VertexPositionZ(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nSharedCluster(0),
  fProton_TPC_nClusterFindable(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_ID(0),
  fProton_Event_Identifier(0),
  fTree_Deuteron(0),
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
  fDeuteron_Event_VertexPositionZ(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nSharedCluster(0),
  fDeuteron_TPC_nClusterFindable(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Identifier(0),
  fTree_AntiProton(0),
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
  fAntiProton_Event_VertexPositionZ(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nSharedCluster(0),
  fAntiProton_TPC_nClusterFindable(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_ID(0),
  fAntiProton_Event_Identifier(0),
  fTree_AntiDeuteron(0),
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
  fAntiDeuteron_Event_VertexPositionZ(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nSharedCluster(0),
  fAntiDeuteron_TPC_nClusterFindable(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Identifier(0)
{


}



AliAnalysisTask_pd::AliAnalysisTask_pd(const char *name,int CollisionSystem) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(CollisionSystem),
  fTree_Proton(0),
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
  fProton_Event_VertexPositionZ(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nSharedCluster(0),
  fProton_TPC_nClusterFindable(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_ID(0),
  fProton_Event_Identifier(0),
  fTree_Deuteron(0),
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
  fDeuteron_Event_VertexPositionZ(0),
  fDeuteron_TPC_nCrossedRows(0),
  fDeuteron_TPC_nSharedCluster(0),
  fDeuteron_TPC_nClusterFindable(0),
  fDeuteron_TPC_nCluster(0),
  fDeuteron_ITS_nCluster(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Identifier(0),
  fTree_AntiProton(0),
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
  fAntiProton_Event_VertexPositionZ(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nSharedCluster(0),
  fAntiProton_TPC_nClusterFindable(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_ID(0),
  fAntiProton_Event_Identifier(0),
  fTree_AntiDeuteron(0),
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
  fAntiDeuteron_Event_VertexPositionZ(0),
  fAntiDeuteron_TPC_nCrossedRows(0),
  fAntiDeuteron_TPC_nSharedCluster(0),
  fAntiDeuteron_TPC_nClusterFindable(0),
  fAntiDeuteron_TPC_nCluster(0),
  fAntiDeuteron_ITS_nCluster(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Identifier(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());

}

  
AliAnalysisTask_pd::~AliAnalysisTask_pd()
{

  if(fTree_Proton)
    {
      delete fTree_Proton;
    }

  if(fTree_Deuteron)
    {
      delete fTree_Deuteron;
    }

  if(fTree_AntiProton)
    {
      delete fTree_AntiProton;
    }

  if(fTree_AntiDeuteron)
    {
      delete fTree_AntiDeuteron;
    }


}







void AliAnalysisTask_pd::UserCreateOutputObjects()
{



  fTree_Proton = new TTree("fTree_Proton","fTree_Proton");
  fTree_Proton->Branch("Proton_px",&fProton_px,"Proton_px/f");
  fTree_Proton->Branch("Proton_py",&fProton_py,"Proton_py/f");
  fTree_Proton->Branch("Proton_pz",&fProton_pz,"Proton_pz/f");
  fTree_Proton->Branch("Proton_pTPC",&fProton_pTPC,"Proton_pTPC/f");
  fTree_Proton->Branch("Proton_Eta",&fProton_Eta,"Proton_Eta/f");
  fTree_Proton->Branch("Proton_Phi",&fProton_Phi,"Proton_Phi/f");
  fTree_Proton->Branch("Proton_TPC_Chi2overNDF",&fProton_TPC_Chi2overNDF,"Proton_TPC_Chi2overNDF/f");
  fTree_Proton->Branch("Proton_TPC_dEdx",&fProton_TPC_dEdx,"Proton_TPC_dEdx/f");
  fTree_Proton->Branch("Proton_TPC_dEdx_Sigma",&fProton_TPC_dEdx_Sigma,"Proton_TPC_dEdx_Sigma/f");
  fTree_Proton->Branch("Proton_TOF_Beta",&fProton_TOF_Beta,"Proton_TOF_Beta/f");
  fTree_Proton->Branch("Proton_TOF_Beta_Sigma",&fProton_TOF_Beta_Sigma,"Proton_TOF_Beta_Sigma/f");
  fTree_Proton->Branch("Proton_TOF_Mass2",&fProton_TOF_Mass2,"Proton_TOF_Mass2/f");
  fTree_Proton->Branch("Proton_TOF_Mass2_Sigma",&fProton_TOF_Mass2_Sigma,"Proton_TOF_Mass2_Sigma/f");
  fTree_Proton->Branch("Proton_ITS_dEdx",&fProton_ITS_dEdx,"Proton_ITS_dEdx/f");
  fTree_Proton->Branch("Proton_ITS_dEdx_Sigma",&fProton_ITS_dEdx_Sigma,"Proton_ITS_dEdx_Sigma/f");
  fTree_Proton->Branch("Proton_DCAxy",&fProton_DCAxy,"Proton_DCAxy/f");
  fTree_Proton->Branch("Proton_DCAz",&fProton_DCAz,"Proton_DCAz/f");
  fTree_Proton->Branch("Proton_Event_Centrality",&fProton_Event_Centrality,"Proton_Event_Centrality/f");
  fTree_Proton->Branch("Proton_Event_VertexPositionZ",&fProton_Event_VertexPositionZ,"Proton_Event_VertexPositionZ/f");
  fTree_Proton->Branch("Proton_TPC_nCrossedRows",&fProton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fTree_Proton->Branch("Proton_TPC_nSharedCluster",&fProton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fTree_Proton->Branch("Proton_TPC_nClusterFindable",&fProton_TPC_nClusterFindable,"Proton_TPC_nClusterFindable/s");
  fTree_Proton->Branch("Proton_TPC_nCluster",&fProton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fTree_Proton->Branch("Proton_ITS_nCluster",&fProton_ITS_nCluster,"Proton_ITS_nCluster/s");
  fTree_Proton->Branch("Proton_ID",&fProton_ID,"Proton_ID/i");
  fTree_Proton->Branch("Proton_Event_Identifier",&fProton_Event_Identifier,"Proton_Event_Identifier/i");


  fTree_Deuteron = new TTree("fTree_Deuteron","fTree_Deuteron");
  fTree_Deuteron->Branch("Deuteron_px",&fDeuteron_px,"Deuteron_px/f");
  fTree_Deuteron->Branch("Deuteron_py",&fDeuteron_py,"Deuteron_py/f");
  fTree_Deuteron->Branch("Deuteron_pz",&fDeuteron_pz,"Deuteron_pz/f");
  fTree_Deuteron->Branch("Deuteron_pTPC",&fDeuteron_pTPC,"Deuteron_pTPC/f");
  fTree_Deuteron->Branch("Deuteron_Eta",&fDeuteron_Eta,"Deuteron_Eta/f");
  fTree_Deuteron->Branch("Deuteron_Phi",&fDeuteron_Phi,"Deuteron_Phi/f");
  fTree_Deuteron->Branch("Deuteron_TPC_Chi2overNDF",&fDeuteron_TPC_Chi2overNDF,"Deuteron_TPC_Chi2overNDF/f");
  fTree_Deuteron->Branch("Deuteron_TPC_dEdx",&fDeuteron_TPC_dEdx,"Deuteron_TPC_dEdx/f");
  fTree_Deuteron->Branch("Deuteron_TPC_dEdx_Sigma",&fDeuteron_TPC_dEdx_Sigma,"Deuteron_TPC_dEdx_Sigma/f");
  fTree_Deuteron->Branch("Deuteron_TOF_Beta",&fDeuteron_TOF_Beta,"Deuteron_TOF_Beta/f");
  fTree_Deuteron->Branch("Deuteron_TOF_Beta_Sigma",&fDeuteron_TOF_Beta_Sigma,"Deuteron_TOF_Beta_Sigma/f");
  fTree_Deuteron->Branch("Deuteron_TOF_Mass2",&fDeuteron_TOF_Mass2,"Deuteron_TOF_Mass2/f");
  fTree_Deuteron->Branch("Deuteron_TOF_Mass2_Sigma",&fDeuteron_TOF_Mass2_Sigma,"Deuteron_TOF_Mass2_Sigma/f");
  fTree_Deuteron->Branch("Deuteron_ITS_dEdx",&fDeuteron_ITS_dEdx,"Deuteron_ITS_dEdx/f");
  fTree_Deuteron->Branch("Deuteron_ITS_dEdx_Sigma",&fDeuteron_ITS_dEdx_Sigma,"Deuteron_ITS_dEdx_Sigma/f");
  fTree_Deuteron->Branch("Deuteron_DCAxy",&fDeuteron_DCAxy,"Deuteron_DCAxy/f");
  fTree_Deuteron->Branch("Deuteron_DCAz",&fDeuteron_DCAz,"Deuteron_DCAz/f");
  fTree_Deuteron->Branch("Deuteron_Event_Centrality",&fDeuteron_Event_Centrality,"Deuteron_Event_Centrality/f");
  fTree_Deuteron->Branch("Deuteron_Event_VertexPositionZ",&fDeuteron_Event_VertexPositionZ,"Deuteron_Event_VertexPositionZ/f");
  fTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&fDeuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&fDeuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fTree_Deuteron->Branch("Deuteron_TPC_nClusterFindable",&fDeuteron_TPC_nClusterFindable,"Deuteron_TPC_nClusterFindable/s");
  fTree_Deuteron->Branch("Deuteron_TPC_nCluster",&fDeuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fTree_Deuteron->Branch("Deuteron_ITS_nCluster",&fDeuteron_ITS_nCluster,"Deuteron_ITS_nCluster/s");
  fTree_Deuteron->Branch("Deuteron_ID",&fDeuteron_ID,"Deuteron_ID/i");
  fTree_Deuteron->Branch("Deuteron_Event_Identifier",&fDeuteron_Event_Identifier,"Deuteron_Event_Identifier/i");


  fTree_AntiProton = new TTree("fTree_AntiProton","fTree_AntiProton");
  fTree_AntiProton->Branch("AntiProton_px",&fAntiProton_px,"AntiProton_px/f");
  fTree_AntiProton->Branch("AntiProton_py",&fAntiProton_py,"AntiProton_py/f");
  fTree_AntiProton->Branch("AntiProton_pz",&fAntiProton_pz,"AntiProton_pz/f");
  fTree_AntiProton->Branch("AntiProton_pTPC",&fAntiProton_pTPC,"AntiProton_pTPC/f");
  fTree_AntiProton->Branch("AntiProton_Eta",&fAntiProton_Eta,"AntiProton_Eta/f");
  fTree_AntiProton->Branch("AntiProton_Phi",&fAntiProton_Phi,"AntiProton_Phi/f");
  fTree_AntiProton->Branch("AntiProton_TPC_Chi2overNDF",&fAntiProton_TPC_Chi2overNDF,"AntiProton_TPC_Chi2overNDF/f");
  fTree_AntiProton->Branch("AntiProton_TPC_dEdx",&fAntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/f");
  fTree_AntiProton->Branch("AntiProton_TPC_dEdx_Sigma",&fAntiProton_TPC_dEdx_Sigma,"AntiProton_TPC_dEdx_Sigma/f");
  fTree_AntiProton->Branch("AntiProton_TOF_Beta",&fAntiProton_TOF_Beta,"AntiProton_TOF_Beta/f");
  fTree_AntiProton->Branch("AntiProton_TOF_Beta_Sigma",&fAntiProton_TOF_Beta_Sigma,"AntiProton_TOF_Beta_Sigma/f");
  fTree_AntiProton->Branch("AntiProton_TOF_Mass2",&fAntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/f");
  fTree_AntiProton->Branch("AntiProton_TOF_Mass2_Sigma",&fAntiProton_TOF_Mass2_Sigma,"AntiProton_TOF_Mass2_Sigma/f");
  fTree_AntiProton->Branch("AntiProton_ITS_dEdx",&fAntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/f");
  fTree_AntiProton->Branch("AntiProton_ITS_dEdx_Sigma",&fAntiProton_ITS_dEdx_Sigma,"AntiProton_ITS_dEdx_Sigma/f");
  fTree_AntiProton->Branch("AntiProton_DCAxy",&fAntiProton_DCAxy,"AntiProton_DCAxy/f");
  fTree_AntiProton->Branch("AntiProton_DCAz",&fAntiProton_DCAz,"AntiProton_DCAz/f");
  fTree_AntiProton->Branch("AntiProton_Event_Centrality",&fAntiProton_Event_Centrality,"AntiProton_Event_Centrality/f");
  fTree_AntiProton->Branch("AntiProton_Event_VertexPositionZ",&fAntiProton_Event_VertexPositionZ,"AntiProton_Event_VertexPositionZ/f");
  fTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&fAntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&fAntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fTree_AntiProton->Branch("AntiProton_TPC_nClusterFindable",&fAntiProton_TPC_nClusterFindable,"AntiProton_TPC_nClusterFindable/s");
  fTree_AntiProton->Branch("AntiProton_TPC_nCluster",&fAntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fTree_AntiProton->Branch("AntiProton_ITS_nCluster",&fAntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/s");
  fTree_AntiProton->Branch("AntiProton_ID",&fAntiProton_ID,"AntiProton_ID/i");
  fTree_AntiProton->Branch("AntiProton_Event_Identifier",&fAntiProton_Event_Identifier,"AntiProton_Event_Identifier/i");


  fTree_AntiDeuteron = new TTree("fTree_AntiDeuteron","fTree_AntiDeuteron");
  fTree_AntiDeuteron->Branch("AntiDeuteron_px",&fAntiDeuteron_px,"AntiDeuteron_px/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_py",&fAntiDeuteron_py,"AntiDeuteron_py/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_pz",&fAntiDeuteron_pz,"AntiDeuteron_pz/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&fAntiDeuteron_pTPC,"AntiDeuteron_pTPC/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&fAntiDeuteron_Eta,"AntiDeuteron_Eta/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&fAntiDeuteron_Phi,"AntiDeuteron_Phi/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2overNDF",&fAntiDeuteron_TPC_Chi2overNDF,"AntiDeuteron_TPC_Chi2overNDF/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&fAntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_Sigma",&fAntiDeuteron_TPC_dEdx_Sigma,"AntiDeuteron_TPC_dEdx_Sigma/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta",&fAntiDeuteron_TOF_Beta,"AntiDeuteron_TOF_Beta/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta_Sigma",&fAntiDeuteron_TOF_Beta_Sigma,"AntiDeuteron_TOF_Beta_Sigma/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&fAntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_Sigma",&fAntiDeuteron_TOF_Mass2_Sigma,"AntiDeuteron_TOF_Mass2_Sigma/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&fAntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_Sigma",&fAntiDeuteron_ITS_dEdx_Sigma,"AntiDeuteron_ITS_dEdx_Sigma/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&fAntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&fAntiDeuteron_DCAz,"AntiDeuteron_DCAz/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_Event_Centrality",&fAntiDeuteron_Event_Centrality,"AntiDeuteron_Event_Centrality/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_Event_VertexPositionZ",&fAntiDeuteron_Event_VertexPositionZ,"AntiDeuteron_Event_VertexPositionZ/f");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&fAntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&fAntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nClusterFindable",&fAntiDeuteron_TPC_nClusterFindable,"AntiDeuteron_TPC_nClusterFindable/s");
  fTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&fAntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&fAntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/s");
  fTree_AntiDeuteron->Branch("AntiDeuteron_ID",&fAntiDeuteron_ID,"AntiDeuteron_ID/i");
  fTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&fAntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/i");








  




  PostData(1,fTree_Proton);
  PostData(2,fTree_Deuteron);
  PostData(3,fTree_AntiProton);
  PostData(4,fTree_AntiDeuteron);



} // end of UserCreateOutputObjects









void AliAnalysisTask_pd::UserExec(Option_t*)
{

  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAODEvent)::Fatal("AliAnalysisTask_pd::UserExec","No AOD event found!");

  fHeader = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  if(!fHeader)::Fatal("AliAnalysisTask_pd::UserExec","No Header found!");

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse)::Fatal("AliAnalysisTask_pd::UserExec","No PIDResponse found!");


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
  if(!PrimaryVertex)::Warning("AliAnalsisTask_pd::UserExec","No AliAODVertex object found!");
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




  // set up particle masses (GeV/c2)
  const double ProtonMass   = 0.9382720;
  const double DeuteronMass = 1.8756129;


  bool ProtonIsSelected	      = false;
  bool DeuteronIsSelected     = false;
  bool AntiProtonIsSelected   = false;
  bool AntiDeuteronIsSelected = false;



  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ proton selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// proton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply proton cuts
    bool PassedProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,true);
    if(!PassedProtonCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    fProton_px			    = Track->Px();
    fProton_py			    = Track->Py();
    fProton_pz			    = Track->Pz();
    fProton_pTPC		    = Track->GetTPCmomentum();
    fProton_Eta			    = Track->Eta();
    fProton_Phi			    = Track->Phi();
    fProton_TPC_Chi2overNDF	    = Track->GetTPCchi2perNDF();
    fProton_TPC_dEdx		    = Track->GetTPCsignal();
    fProton_TPC_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    fProton_TOF_Beta		    = Track->GetTOFsignal();
    fProton_TOF_Beta_Sigma	    = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
    fProton_TOF_Mass2		    = CalculateMassSquareTOF(*Track);
    fProton_ITS_dEdx		    = Track->GetITSsignal();
    fProton_ITS_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kProton);
    fProton_DCAxy		    = DCAxy;
    fProton_DCAz		    = DCAz;
    fProton_Event_Centrality	    = Centrality;
    fProton_Event_VertexPositionZ   = PrimaryVertexZ;
    fProton_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    fProton_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    fProton_TPC_nClusterFindable    = Track->GetTPCNclsF();
    fProton_TPC_nCluster	    = Track->GetTPCNcls();
    fProton_ITS_nCluster	    = Track->GetITSNcls();
    fProton_ID			    = Track->GetID();
    fProton_Event_Identifier	    = BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    ProtonIsSelected = true;

  } // end of proton loop




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ deuteron selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// deuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply deuteron cuts
    bool PassedDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,true);
    if(!PassedDeuteronCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    fDeuteron_px		    = Track->Px();
    fDeuteron_py		    = Track->Py();
    fDeuteron_pz		    = Track->Pz();
    fDeuteron_pTPC		    = Track->GetTPCmomentum();
    fDeuteron_Eta		    = Track->Eta();
    fDeuteron_Phi		    = Track->Phi();
    fDeuteron_TPC_Chi2overNDF	    = Track->GetTPCchi2perNDF();
    fDeuteron_TPC_dEdx		    = Track->GetTPCsignal();
    fDeuteron_TPC_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    fDeuteron_TOF_Beta		    = Track->GetTOFsignal();
    fDeuteron_TOF_Beta_Sigma	    = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
    fDeuteron_TOF_Mass2		    = CalculateMassSquareTOF(*Track);
    fDeuteron_TOF_Mass2_Sigma	    = CalculateDeuteronSigmaMassSquareTOF(Track->Pt(),fDeuteron_TOF_Mass2,true);
    fDeuteron_ITS_dEdx		    = Track->GetITSsignal();
    fDeuteron_ITS_dEdx_Sigma	    = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kDeuteron);
    fDeuteron_DCAxy		    = DCAxy;
    fDeuteron_DCAz		    = DCAz;
    fDeuteron_Event_Centrality	    = Centrality;
    fDeuteron_Event_VertexPositionZ = PrimaryVertexZ;
    fDeuteron_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    fDeuteron_TPC_nSharedCluster    = Track->GetTPCnclsS();
    fDeuteron_TPC_nClusterFindable  = Track->GetTPCNclsF();
    fDeuteron_TPC_nCluster	    = Track->GetTPCNcls();
    fDeuteron_ITS_nCluster	    = Track->GetITSNcls();
    fDeuteron_ID		    = Track->GetID();
    fDeuteron_Event_Identifier	    = BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);

 
    DeuteronIsSelected = true;

  } // end of deuteron loop


  if(ProtonIsSelected && DeuteronIsSelected)
  {

    fTree_Proton->Fill();
    fTree_Deuteron->Fill();

  } // end of ProtonIsSelected && DeuteronIsSelected








  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antiproton selection loop +++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// antiproton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antiproton cuts
    bool PassedAntiProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,false);
    if(!PassedAntiProtonCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    fAntiProton_px		      = Track->Px();
    fAntiProton_py		      = Track->Py();
    fAntiProton_pz		      = Track->Pz();
    fAntiProton_pTPC		      = Track->GetTPCmomentum();
    fAntiProton_Eta		      = Track->Eta();
    fAntiProton_Phi		      = Track->Phi();
    fAntiProton_TPC_Chi2overNDF	      = Track->GetTPCchi2perNDF();
    fAntiProton_TPC_dEdx	      = Track->GetTPCsignal();
    fAntiProton_TPC_dEdx_Sigma	      = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    fAntiProton_TOF_Beta	      = Track->GetTOFsignal();
    fAntiProton_TOF_Beta_Sigma	      = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
    fAntiProton_TOF_Mass2	      = CalculateMassSquareTOF(*Track);
    fAntiProton_ITS_dEdx	      = Track->GetITSsignal();
    fAntiProton_ITS_dEdx_Sigma	      = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kProton);
    fAntiProton_DCAxy		      = DCAxy;
    fAntiProton_DCAz		      = DCAz;
    fAntiProton_Event_Centrality      = Centrality;
    fAntiProton_Event_VertexPositionZ = PrimaryVertexZ;
    fAntiProton_TPC_nCrossedRows      = Track->GetTPCCrossedRows();
    fAntiProton_TPC_nSharedCluster    = Track->GetTPCnclsS();
    fAntiProton_TPC_nClusterFindable  = Track->GetTPCNclsF();
    fAntiProton_TPC_nCluster	      = Track->GetTPCNcls();
    fAntiProton_ITS_nCluster	      = Track->GetITSNcls();
    fAntiProton_ID		      = Track->GetID();
    fAntiProton_Event_Identifier      = BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    AntiProtonIsSelected = true;

  } // end of antiproton loop




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antideuteron selection loop +++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// antideuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antideuteron cuts
    bool PassedAntiDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,false);
    if(!PassedAntiDeuteronCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    fAntiDeuteron_px			= Track->Px();
    fAntiDeuteron_py			= Track->Py();
    fAntiDeuteron_pz			= Track->Pz();
    fAntiDeuteron_pTPC			= Track->GetTPCmomentum();
    fAntiDeuteron_Eta			= Track->Eta();
    fAntiDeuteron_Phi			= Track->Phi();
    fAntiDeuteron_TPC_Chi2overNDF	= Track->GetTPCchi2perNDF();
    fAntiDeuteron_TPC_dEdx		= Track->GetTPCsignal();
    fAntiDeuteron_TPC_dEdx_Sigma	= fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    fAntiDeuteron_TOF_Beta		= Track->GetTOFsignal();
    fAntiDeuteron_TOF_Beta_Sigma	= fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
    fAntiDeuteron_TOF_Mass2		= CalculateMassSquareTOF(*Track);
    fAntiDeuteron_TOF_Mass2_Sigma	= CalculateDeuteronSigmaMassSquareTOF(Track->Pt(),fAntiDeuteron_TOF_Mass2,false);
    fAntiDeuteron_ITS_dEdx		= Track->GetITSsignal();
    fAntiDeuteron_ITS_dEdx_Sigma	= fPIDResponse->NumberOfSigmasITS(Track,AliPID::kDeuteron);
    fAntiDeuteron_DCAxy			= DCAxy;
    fAntiDeuteron_DCAz			= DCAz;
    fAntiDeuteron_Event_Centrality	= Centrality;
    fAntiDeuteron_Event_VertexPositionZ	= PrimaryVertexZ;
    fAntiDeuteron_TPC_nCrossedRows	= Track->GetTPCCrossedRows();
    fAntiDeuteron_TPC_nSharedCluster    = Track->GetTPCnclsS();
    fAntiDeuteron_TPC_nClusterFindable  = Track->GetTPCNclsF();
    fAntiDeuteron_TPC_nCluster		= Track->GetTPCNcls();
    fAntiDeuteron_ITS_nCluster		= Track->GetITSNcls();
    fAntiDeuteron_ID			= Track->GetID();
    fAntiDeuteron_Event_Identifier	= BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    AntiDeuteronIsSelected = true;

  } // end of antideuteron loop


  if(AntiProtonIsSelected && AntiDeuteronIsSelected)
  {

    fTree_AntiProton->Fill();
    fTree_AntiDeuteron->Fill();

  } // end of AntiProtonIsSelected && AntiDeuteronIsSelected









  PostData(1,fTree_Proton);
  PostData(2,fTree_Deuteron);
  PostData(3,fTree_AntiProton);
  PostData(4,fTree_AntiDeuteron);

} // end of UserExec















void AliAnalysisTask_pd::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
double AliAnalysisTask_pd::CalculateBetaTOF(AliAODTrack &track)
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
double AliAnalysisTask_pd::CalculateMassSquareTOF(AliAODTrack &track)
{

  double p = track.P();
  double beta = CalculateBetaTOF(track);
  double mass2 = -999.0;

  if(beta > 0.0){

    mass2 = (1/(beta*beta)-1) * (p*p);

  }

  return mass2;

} // end of CalculateMassSquareTOF








double AliAnalysisTask_pd::CalculateDeuteronSigmaMassSquareTOF(double pT, double massSq, bool isMatter)
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
bool AliAnalysisTask_pd::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter)
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
bool AliAnalysisTask_pd::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter)
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
  double Deuteron_ITS_nSigma_max = 60;


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
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse.CheckPIDStatus(AliPIDResponse::kITS,&Track);
  bool ITSisOK = false;
  if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;

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



