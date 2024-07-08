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

#include "AliAnalysisTask_pd_CreateTrees_PairsOnly.h"

using namespace std;
ClassImp(AliAnalysisTask_pd_CreateTrees_PairsOnly) 






AliAnalysisTask_pd_CreateTrees_PairsOnly::AliAnalysisTask_pd_CreateTrees_PairsOnly() : AliAnalysisTaskSE(),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(0),
  fUseOpenCuts(0),
  fIsMC(0),
  fSavePairsOnly(0),
  fTimeRangeCut(),
  fSaveTree_Proton(0),
  fProton_px(-999.0),
  fProton_py(-999.0),
  fProton_pz(-999.0),
  fProton_px_Generated(-999.0),
  fProton_py_Generated(-999.0),
  fProton_pz_Generated(-999.0),
  fProton_pTPC(-999.0),
  fProton_Eta(-999.0),
  fProton_Phi(-999.0),
  fProton_TPC_Chi2(-999.0),
  fProton_TPC_dEdx(-999.0),
  fProton_TPC_dEdx_Sigma(-999.0),
  fProton_TOF_Mass2(-999.0),
  fProton_TOF_Mass2_Sigma(-999.0),
  fProton_ITS_dEdx(-999.0),
  fProton_ITS_dEdx_Sigma(-999.0),
  fProton_DCAxy(-999.0),
  fProton_DCAz(-999.0),
  fProton_Event_Centrality(-999.0),
  fProton_Event_PrimaryVertexZ(-999.0),
  fProton_Event_BField(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nFindableCluster(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_PDG(0),
  fProton_MotherPDG(0),
  fProton_FilterBit(0),
  fProton_Event_Multiplicity(0),
  fProton_Event_Identifier(0),
  fProton_Event_RunNumber(0),
  fProton_ITS_Layer0(0),
  fProton_ITS_Layer1(0),
  fProton_ITS_Layer2(0),
  fProton_ITS_Layer3(0),
  fProton_ITS_Layer4(0),
  fProton_ITS_Layer5(0),
  fProton_Event_IsFirstParticle(0),
  fProton_Event_TimeStamp(0),
  fProton_Event_RandomCrossCheckNumber(0),
  fProton_ID(0),
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
  fDeuteron_Event_TimeStamp(0),
  fDeuteron_Event_RandomCrossCheckNumber(0),
  fDeuteron_ID(0),
  fSaveTree_AntiProton(0),
  fAntiProton_px(-999.0),
  fAntiProton_py(-999.0),
  fAntiProton_pz(-999.0),
  fAntiProton_px_Generated(-999.0),
  fAntiProton_py_Generated(-999.0),
  fAntiProton_pz_Generated(-999.0),
  fAntiProton_pTPC(-999.0),
  fAntiProton_Eta(-999.0),
  fAntiProton_Phi(-999.0),
  fAntiProton_TPC_Chi2(-999.0),
  fAntiProton_TPC_dEdx(-999.0),
  fAntiProton_TPC_dEdx_Sigma(-999.0),
  fAntiProton_TOF_Mass2(-999.0),
  fAntiProton_TOF_Mass2_Sigma(-999.0),
  fAntiProton_ITS_dEdx(-999.0),
  fAntiProton_ITS_dEdx_Sigma(-999.0),
  fAntiProton_DCAxy(-999.0),
  fAntiProton_DCAz(-999.0),
  fAntiProton_Event_Centrality(-999.0),
  fAntiProton_Event_PrimaryVertexZ(-999.0),
  fAntiProton_Event_BField(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nFindableCluster(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_PDG(0),
  fAntiProton_MotherPDG(0),
  fAntiProton_FilterBit(0),
  fAntiProton_Event_Multiplicity(0),
  fAntiProton_Event_Identifier(0),
  fAntiProton_Event_RunNumber(0),
  fAntiProton_ITS_Layer0(0),
  fAntiProton_ITS_Layer1(0),
  fAntiProton_ITS_Layer2(0),
  fAntiProton_ITS_Layer3(0),
  fAntiProton_ITS_Layer4(0),
  fAntiProton_ITS_Layer5(0),
  fAntiProton_Event_IsFirstParticle(0),
  fAntiProton_Event_TimeStamp(0),
  fAntiProton_Event_RandomCrossCheckNumber(0),
  fAntiProton_ID(0),
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
  fAntiDeuteron_Event_IsFirstParticle(0),
  fAntiDeuteron_Event_TimeStamp(0),
  fAntiDeuteron_Event_RandomCrossCheckNumber(0),
  fAntiDeuteron_ID(0)
{


}



AliAnalysisTask_pd_CreateTrees_PairsOnly::AliAnalysisTask_pd_CreateTrees_PairsOnly(const char *name,Int_t CollisionSystem, Bool_t UseOpenCuts, Bool_t isMC, Bool_t SavePairsOnly) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(CollisionSystem),
  fUseOpenCuts(UseOpenCuts),
  fIsMC(isMC),
  fSavePairsOnly(SavePairsOnly),
  fTimeRangeCut(),
  fSaveTree_Proton(0),
  fProton_px(-999.0),
  fProton_py(-999.0),
  fProton_pz(-999.0),
  fProton_px_Generated(-999.0),
  fProton_py_Generated(-999.0),
  fProton_pz_Generated(-999.0),
  fProton_pTPC(-999.0),
  fProton_Eta(-999.0),
  fProton_Phi(-999.0),
  fProton_TPC_Chi2(-999.0),
  fProton_TPC_dEdx(-999.0),
  fProton_TPC_dEdx_Sigma(-999.0),
  fProton_TOF_Mass2(-999.0),
  fProton_TOF_Mass2_Sigma(-999.0),
  fProton_ITS_dEdx(-999.0),
  fProton_ITS_dEdx_Sigma(-999.0),
  fProton_DCAxy(-999.0),
  fProton_DCAz(-999.0),
  fProton_Event_Centrality(-999.0),
  fProton_Event_PrimaryVertexZ(-999.0),
  fProton_Event_BField(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nFindableCluster(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_PDG(0),
  fProton_MotherPDG(0),
  fProton_FilterBit(0),
  fProton_Event_Multiplicity(0),
  fProton_Event_Identifier(0),
  fProton_Event_RunNumber(0),
  fProton_ITS_Layer0(0),
  fProton_ITS_Layer1(0),
  fProton_ITS_Layer2(0),
  fProton_ITS_Layer3(0),
  fProton_ITS_Layer4(0),
  fProton_ITS_Layer5(0),
  fProton_Event_IsFirstParticle(0),
  fProton_Event_TimeStamp(0),
  fProton_Event_RandomCrossCheckNumber(0),
  fProton_ID(0),
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
  fDeuteron_Event_TimeStamp(0),
  fDeuteron_Event_RandomCrossCheckNumber(0),
  fDeuteron_ID(0),
  fSaveTree_AntiProton(0),
  fAntiProton_px(-999.0),
  fAntiProton_py(-999.0),
  fAntiProton_pz(-999.0),
  fAntiProton_px_Generated(-999.0),
  fAntiProton_py_Generated(-999.0),
  fAntiProton_pz_Generated(-999.0),
  fAntiProton_pTPC(-999.0),
  fAntiProton_Eta(-999.0),
  fAntiProton_Phi(-999.0),
  fAntiProton_TPC_Chi2(-999.0),
  fAntiProton_TPC_dEdx(-999.0),
  fAntiProton_TPC_dEdx_Sigma(-999.0),
  fAntiProton_TOF_Mass2(-999.0),
  fAntiProton_TOF_Mass2_Sigma(-999.0),
  fAntiProton_ITS_dEdx(-999.0),
  fAntiProton_ITS_dEdx_Sigma(-999.0),
  fAntiProton_DCAxy(-999.0),
  fAntiProton_DCAz(-999.0),
  fAntiProton_Event_Centrality(-999.0),
  fAntiProton_Event_PrimaryVertexZ(-999.0),
  fAntiProton_Event_BField(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nFindableCluster(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_PDG(0),
  fAntiProton_MotherPDG(0),
  fAntiProton_FilterBit(0),
  fAntiProton_Event_Multiplicity(0),
  fAntiProton_Event_Identifier(0),
  fAntiProton_Event_RunNumber(0),
  fAntiProton_ITS_Layer0(0),
  fAntiProton_ITS_Layer1(0),
  fAntiProton_ITS_Layer2(0),
  fAntiProton_ITS_Layer3(0),
  fAntiProton_ITS_Layer4(0),
  fAntiProton_ITS_Layer5(0),
  fAntiProton_Event_IsFirstParticle(0),
  fAntiProton_Event_TimeStamp(0),
  fAntiProton_Event_RandomCrossCheckNumber(0),
  fAntiProton_ID(0),
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
  fAntiDeuteron_Event_IsFirstParticle(0),
  fAntiDeuteron_Event_TimeStamp(0),
  fAntiDeuteron_Event_RandomCrossCheckNumber(0),
  fAntiDeuteron_ID(0)
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
  fSaveTree_Proton->Branch("Proton_px",&fProton_px,"Proton_px/F");
  fSaveTree_Proton->Branch("Proton_py",&fProton_py,"Proton_py/F");
  fSaveTree_Proton->Branch("Proton_pz",&fProton_pz,"Proton_pz/F");
  if(fIsMC == true){
  fSaveTree_Proton->Branch("Proton_px_Generated",&fProton_px_Generated,"Proton_px_Generated/F");
  fSaveTree_Proton->Branch("Proton_py_Generated",&fProton_py_Generated,"Proton_py_Generated/F");
  fSaveTree_Proton->Branch("Proton_pz_Generated",&fProton_pz_Generated,"Proton_pz_Generated/F");
  }
  fSaveTree_Proton->Branch("Proton_pTPC",&fProton_pTPC,"Proton_pTPC/F");
  fSaveTree_Proton->Branch("Proton_Eta",&fProton_Eta,"Proton_Eta/F");
  fSaveTree_Proton->Branch("Proton_Phi",&fProton_Phi,"Proton_Phi/F");
  fSaveTree_Proton->Branch("Proton_TPC_Chi2",&fProton_TPC_Chi2,"Proton_TPC_Chi2/F");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx",&fProton_TPC_dEdx,"Proton_TPC_dEdx/F");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx_Sigma",&fProton_TPC_dEdx_Sigma,"Proton_TPC_dEdx_Sigma/F");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2",&fProton_TOF_Mass2,"Proton_TOF_Mass2/F");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2_Sigma",&fProton_TOF_Mass2_Sigma,"Proton_TOF_Mass2_Sigma/F");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx",&fProton_ITS_dEdx,"Proton_ITS_dEdx/F");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx_Sigma",&fProton_ITS_dEdx_Sigma,"Proton_ITS_dEdx_Sigma/F");
  fSaveTree_Proton->Branch("Proton_DCAxy",&fProton_DCAxy,"Proton_DCAxy/F");
  fSaveTree_Proton->Branch("Proton_DCAz",&fProton_DCAz,"Proton_DCAz/F");
  fSaveTree_Proton->Branch("Proton_Event_Centrality",&fProton_Event_Centrality,"Proton_Event_Centrality/F");
  fSaveTree_Proton->Branch("Proton_Event_PrimaryVertexZ",&fProton_Event_PrimaryVertexZ,"Proton_Event_PrimaryVertexZ/F");
  fSaveTree_Proton->Branch("Proton_Event_BField",&fProton_Event_BField,"Proton_Event_BField/O");
  fSaveTree_Proton->Branch("Proton_TPC_nCrossedRows",&fProton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fSaveTree_Proton->Branch("Proton_TPC_nFindableCluster",&fProton_TPC_nFindableCluster,"Proton_TPC_nFindableCluster/s");
  fSaveTree_Proton->Branch("Proton_TPC_nCluster",&fProton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fSaveTree_Proton->Branch("Proton_ITS_nCluster",&fProton_ITS_nCluster,"Proton_ITS_nCluster/s");
  if(fIsMC == true){
  fSaveTree_Proton->Branch("Proton_PDG",&fProton_PDG,"Proton_PDG/I");
  fSaveTree_Proton->Branch("Proton_MotherPDG",&fProton_MotherPDG,"Proton_MotherPDG/I");
  }
  fSaveTree_Proton->Branch("Proton_FilterBit",&fProton_FilterBit,"Proton_FilterBit/O");
  fSaveTree_Proton->Branch("Proton_Event_Multiplicity",&fProton_Event_Multiplicity,"Proton_Event_Multiplicity/I");
  fSaveTree_Proton->Branch("Proton_Event_Identifier",&fProton_Event_Identifier,"Proton_Event_Identifier/l");
  fSaveTree_Proton->Branch("Proton_Event_RunNumber",&fProton_Event_RunNumber,"Proton_Event_RunNumber/I");
  fSaveTree_Proton->Branch("Proton_ITS_Layer0",&fProton_ITS_Layer0,"Proton_ITS_Layer0/O");
  fSaveTree_Proton->Branch("Proton_ITS_Layer1",&fProton_ITS_Layer1,"Proton_ITS_Layer1/O");
  fSaveTree_Proton->Branch("Proton_ITS_Layer2",&fProton_ITS_Layer2,"Proton_ITS_Layer2/O");
  fSaveTree_Proton->Branch("Proton_ITS_Layer3",&fProton_ITS_Layer3,"Proton_ITS_Layer3/O");
  fSaveTree_Proton->Branch("Proton_ITS_Layer4",&fProton_ITS_Layer4,"Proton_ITS_Layer4/O");
  fSaveTree_Proton->Branch("Proton_ITS_Layer5",&fProton_ITS_Layer5,"Proton_ITS_Layer5/O");
  fSaveTree_Proton->Branch("Proton_Event_IsFirstParticle",&fProton_Event_IsFirstParticle,"Proton_Event_IsFirstParticle/O");
  fSaveTree_Proton->Branch("Proton_Event_TimeStamp",&fProton_Event_TimeStamp,"Proton_Event_TimeStamp/i");
  fSaveTree_Proton->Branch("Proton_Event_RandomCrossCheckNumber",&fProton_Event_RandomCrossCheckNumber,"Proton_Event_RandomCrossCheckNumber/i");
  fSaveTree_Proton->Branch("Proton_ID",&fProton_ID,"Proton_ID/I");




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
  fSaveTree_Deuteron->Branch("Deuteron_Event_TimeStamp",&fDeuteron_Event_TimeStamp,"Deuteron_Event_TimeStamp/i");
  fSaveTree_Deuteron->Branch("Deuteron_Event_RandomCrossCheckNumber",&fDeuteron_Event_RandomCrossCheckNumber,"Deuteron_Event_RandomCrossCheckNumber/i");
  fSaveTree_Deuteron->Branch("Deuteron_ID",&fDeuteron_ID,"Deuteron_ID/I");



  fSaveTree_AntiProton = new TTree("fSaveTree_AntiProton","fSaveTree_AntiProton");
  fSaveTree_AntiProton->Branch("AntiProton_px",&fAntiProton_px,"AntiProton_px/F");
  fSaveTree_AntiProton->Branch("AntiProton_py",&fAntiProton_py,"AntiProton_py/F");
  fSaveTree_AntiProton->Branch("AntiProton_pz",&fAntiProton_pz,"AntiProton_pz/F");
  if(fIsMC == true){
  fSaveTree_AntiProton->Branch("AntiProton_px_Generated",&fAntiProton_px_Generated,"AntiProton_px_Generated/F");
  fSaveTree_AntiProton->Branch("AntiProton_py_Generated",&fAntiProton_py_Generated,"AntiProton_py_Generated/F");
  fSaveTree_AntiProton->Branch("AntiProton_pz_Generated",&fAntiProton_pz_Generated,"AntiProton_pz_Generated/F");
  }
  fSaveTree_AntiProton->Branch("AntiProton_pTPC",&fAntiProton_pTPC,"AntiProton_pTPC/F");
  fSaveTree_AntiProton->Branch("AntiProton_Eta",&fAntiProton_Eta,"AntiProton_Eta/F");
  fSaveTree_AntiProton->Branch("AntiProton_Phi",&fAntiProton_Phi,"AntiProton_Phi/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_Chi2",&fAntiProton_TPC_Chi2,"AntiProton_TPC_Chi2/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx",&fAntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx_Sigma",&fAntiProton_TPC_dEdx_Sigma,"AntiProton_TPC_dEdx_Sigma/F");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2",&fAntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/F");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2_Sigma",&fAntiProton_TOF_Mass2_Sigma,"AntiProton_TOF_Mass2_Sigma/F");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx",&fAntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/F");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx_Sigma",&fAntiProton_ITS_dEdx_Sigma,"AntiProton_ITS_dEdx_Sigma/F");
  fSaveTree_AntiProton->Branch("AntiProton_DCAxy",&fAntiProton_DCAxy,"AntiProton_DCAxy/F");
  fSaveTree_AntiProton->Branch("AntiProton_DCAz",&fAntiProton_DCAz,"AntiProton_DCAz/F");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Centrality",&fAntiProton_Event_Centrality,"AntiProton_Event_Centrality/F");
  fSaveTree_AntiProton->Branch("AntiProton_Event_PrimaryVertexZ",&fAntiProton_Event_PrimaryVertexZ,"AntiProton_Event_PrimaryVertexZ/F");
  fSaveTree_AntiProton->Branch("AntiProton_Event_BField",&fAntiProton_Event_BField,"AntiProton_Event_BField/O");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&fAntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nFindableCluster",&fAntiProton_TPC_nFindableCluster,"AntiProton_TPC_nFindableCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCluster",&fAntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_nCluster",&fAntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/s");
  if(fIsMC == true){
  fSaveTree_AntiProton->Branch("AntiProton_PDG",&fAntiProton_PDG,"AntiProton_PDG/I");
  fSaveTree_AntiProton->Branch("AntiProton_MotherPDG",&fAntiProton_MotherPDG,"AntiProton_MotherPDG/I");
  }
  fSaveTree_AntiProton->Branch("AntiProton_FilterBit",&fAntiProton_FilterBit,"AntiProton_FilterBit/O");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Multiplicity",&fAntiProton_Event_Multiplicity,"AntiProton_Event_Multiplicity/I");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Identifier",&fAntiProton_Event_Identifier,"AntiProton_Event_Identifier/l");
  fSaveTree_AntiProton->Branch("AntiProton_Event_RunNumber",&fAntiProton_Event_RunNumber,"AntiProton_Event_RunNumber/I");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_Layer0",&fAntiProton_ITS_Layer0,"AntiProton_ITS_Layer0/O");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_Layer1",&fAntiProton_ITS_Layer1,"AntiProton_ITS_Layer1/O");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_Layer2",&fAntiProton_ITS_Layer2,"AntiProton_ITS_Layer2/O");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_Layer3",&fAntiProton_ITS_Layer3,"AntiProton_ITS_Layer3/O");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_Layer4",&fAntiProton_ITS_Layer4,"AntiProton_ITS_Layer4/O");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_Layer5",&fAntiProton_ITS_Layer5,"AntiProton_ITS_Layer5/O");
  fSaveTree_AntiProton->Branch("AntiProton_Event_IsFirstParticle",&fAntiProton_Event_IsFirstParticle,"AntiProton_Event_IsFirstParticle/O");
  fSaveTree_AntiProton->Branch("AntiProton_Event_TimeStamp",&fAntiProton_Event_TimeStamp,"AntiProton_Event_TimeStamp/i");
  fSaveTree_AntiProton->Branch("AntiProton_Event_RandomCrossCheckNumber",&fAntiProton_Event_RandomCrossCheckNumber,"AntiProton_Event_RandomCrossCheckNumber/i");
  fSaveTree_AntiProton->Branch("AntiProton_ID",&fAntiProton_ID,"AntiProton_ID/I");



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
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_TimeStamp",&fAntiDeuteron_Event_TimeStamp,"AntiDeuteron_Event_TimeStamp/i");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_RandomCrossCheckNumber",&fAntiDeuteron_Event_RandomCrossCheckNumber,"AntiDeuteron_Event_RandomCrossCheckNumber/i");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ID",&fAntiDeuteron_ID,"AntiDeuteron_ID/I");




  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiProton);
  PostData(4,fSaveTree_AntiDeuteron);



} // end of UserCreateOutputObjects









void AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec(Option_t*)
{

//  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAODEvent)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AOD event found!");

  fHeader = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  if(!fHeader)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No Header found!");

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No PIDResponse found!");

  if(fIsMC == true){
  
    fMCEvent = MCEvent(); 
    if(!fMCEvent)::Fatal("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No MC event found!");

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
  if(!PrimaryVertex)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnlyd::UserExec","No AliAODVertex object found!");
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

  const char *CurrentFileName = AliAnalysisTaskSE::CurrentFileName();


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

  // get random number between 0 and 65535 (maximum of unsigned short)
  UInt_t RandomCrossCheckNumber = (UInt_t) RandomGenerator->Integer(4294967295);

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
    cout << "CurrentFileName:\t\t" << CurrentFileName << std::endl;
    cout << "PeriodNumber:\t\t\t" << PeriodNumber << endl;
    cout << "TimeStamp:\t\t\t" << TimeStamp << endl;
    cout << "RunNumber:\t\t\t" << RunNumber << endl;
    cout << "OrbitNumber:\t\t\t" << OrbitNumber << endl;
    cout << "BunchCrossNumber:\t\t" << BunchCrossNumber << endl;
    cout << "Unique Event ID:\t\t" << EventID << endl;
    cout << "RandomCrossCheckNumber:\t\t" << RandomCrossCheckNumber << std::endl;
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
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");

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
  // +++ Proton selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Int_t nProtonsSelected = 0;
  Int_t nAntiProtonsSelected = 0;
  std::vector<Int_t> ProtonVector;

  // Proton / AntiProton loop
  for(Int_t track = 0; track < nTracks; track++){

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");

    Short_t Charge = 0;
    Charge = Track->Charge();
    if(!(TMath::Abs(Charge) == 1)) continue;

    Bool_t ProtonIsMatter = false;
    if(Charge == +1) ProtonIsMatter = true;
    if(Charge == -1) ProtonIsMatter = false;

    Int_t ParticleSpecies = 0;
    if(ProtonIsMatter == true) ParticleSpecies = 1;
    if(ProtonIsMatter == false) ParticleSpecies = 3;

    Bool_t PassedProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,ParticleSpecies,RunNumber);
    if(PassedProtonCuts == false) continue;

    if(ProtonIsMatter == true){

      ProtonVector.push_back(track);
      nProtonsSelected++;

    } // end of isProton

    if(ProtonIsMatter == false){

      ProtonVector.push_back(track);
      nAntiProtonsSelected++;

    } // end of is AntiProton


  } // end of Proton / AntiProton selection loop


  // if there are no Protons AND no AntiProtons get rid of this event
  if((fSavePairsOnly == true) && (nProtonsSelected == 0) && (nAntiProtonsSelected == 0)) return;


















  Bool_t SaveMatter = false;
  if((fSavePairsOnly == true) && (nDeuteronsSelected > 0) && (nProtonsSelected > 0)) SaveMatter = true;
  if((fSavePairsOnly == false) && ((nDeuteronsSelected > 0) || (nProtonsSelected > 0))) SaveMatter = true;

  Bool_t SaveAntiMatter = false;
  if((fSavePairsOnly == true) && (nAntiDeuteronsSelected > 0) && (nAntiProtonsSelected > 0)) SaveAntiMatter = true;
  if((fSavePairsOnly == false) && ((nAntiDeuteronsSelected > 0) || (nAntiProtonsSelected > 0))) SaveAntiMatter = true;


  if((SaveMatter == true) || (SaveAntiMatter == true)){


    Int_t nDeuteronsInVector = DeuteronVector.size();
    Bool_t IsFirstDeuteron = true;
    Bool_t IsFirstAntideuteron = true;

    // loop over Deuterons / Antideuterons
    for(Int_t iDeuteron = 0; iDeuteron < nDeuteronsInVector; iDeuteron++){

      Int_t track = DeuteronVector.at(iDeuteron);

      AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
      if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");


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
	if(fIsMC == false) fDeuteron_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
	if(fIsMC == true) fDeuteron_TPC_dEdx_Sigma = CalculateSigmadEdxTPC(*Track,2);
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
	fDeuteron_Event_TimeStamp = TimeStamp;
	fDeuteron_Event_RandomCrossCheckNumber = RandomCrossCheckNumber;
	fDeuteron_ID = track;

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
	if(fIsMC == false) fAntiDeuteron_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
	if(fIsMC == true) fAntiDeuteron_TPC_dEdx_Sigma = CalculateSigmadEdxTPC(*Track,4);
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
	fAntiDeuteron_Event_TimeStamp = TimeStamp;
	fAntiDeuteron_Event_RandomCrossCheckNumber = RandomCrossCheckNumber;
	fAntiDeuteron_ID = track;

	fSaveTree_AntiDeuteron->Fill();

	IsFirstAntideuteron = false;

      } // end of isAntideuteron


    } // end of loop over Deuterons / Antideuterons






    Int_t nProtonsInVector = ProtonVector.size();
    Bool_t IsFirstProton = true;
    Bool_t IsFirstAntiProton = true;

    // loop over Protons / Antideuterons
    for(Int_t iProton = 0; iProton < nProtonsInVector; iProton++){

      Int_t track = ProtonVector.at(iProton);

      AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
      if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");


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

      Bool_t ProtonIsMatter = false;
      if(Charge > 0) ProtonIsMatter = true;
      if(Charge < 0) ProtonIsMatter = false;

      Int_t ParticleSpecies = 0;
      if(ProtonIsMatter == true) ParticleSpecies = 1;
      if(ProtonIsMatter == false) ParticleSpecies = 3;


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


      if((ProtonIsMatter == true) && (SaveMatter == true)){

	fProton_px = Track->Px();
	fProton_py = Track->Py();
	fProton_pz = Track->Pz();
	fProton_px_Generated = Generated_px;
	fProton_py_Generated = Generated_py;
	fProton_pz_Generated = Generated_pz;
	fProton_pTPC = Track->GetTPCmomentum();
	fProton_Eta = Track->Eta();
	fProton_Phi = Track->Phi();
	fProton_TPC_Chi2 = Track->GetTPCchi2();
	fProton_TPC_dEdx = Track->GetTPCsignal();
	if(fIsMC == false) fProton_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
	if(fIsMC == true) fProton_TPC_dEdx_Sigma = CalculateSigmadEdxTPC(*Track,1);
	fProton_TOF_Mass2 = (Float_t)TOF_m2;
	fProton_TOF_Mass2_Sigma = (Float_t)TOF_m2_Sigma;
	fProton_ITS_dEdx = (Float_t)ITS_dEdx;
	fProton_ITS_dEdx_Sigma = (Float_t)ITS_dEdx_Sigma;
	fProton_DCAxy = DCAxy;
	fProton_DCAz = DCAz;
	fProton_TPC_nCrossedRows = Track->GetTPCCrossedRows();
	fProton_TPC_nFindableCluster = Track->GetTPCNclsF();
	fProton_TPC_nCluster = Track->GetTPCNcls();
	fProton_ITS_nCluster = (UShort_t)ITS_nCluster;
	fProton_PDG = PDG;
	fProton_MotherPDG = MotherPDG;
	fProton_FilterBit = Track->TestFilterBit(BIT(0));
	fProton_ITS_Layer0 = Track->HasPointOnITSLayer(0);
	fProton_ITS_Layer1 = Track->HasPointOnITSLayer(1);
	fProton_ITS_Layer2 = Track->HasPointOnITSLayer(2);
	fProton_ITS_Layer3 = Track->HasPointOnITSLayer(3);
	fProton_ITS_Layer4 = Track->HasPointOnITSLayer(4);
	fProton_ITS_Layer5 = Track->HasPointOnITSLayer(5);
	fProton_Event_Multiplicity = Multiplicity;
	fProton_Event_Centrality = Centrality;
	fProton_Event_PrimaryVertexZ = PrimaryVertexZ;
	fProton_Event_BField = IsPositiveBFieldPolarity;
	fProton_Event_Identifier = EventID;
	fProton_Event_RunNumber = RunNumber;
	fProton_Event_IsFirstParticle = IsFirstProton;
	fProton_Event_TimeStamp = TimeStamp;
	fProton_Event_RandomCrossCheckNumber = RandomCrossCheckNumber;
	fProton_ID = track;

	fSaveTree_Proton->Fill();

	IsFirstProton = false;

      } // end of isProton



      if((ProtonIsMatter == false) && (SaveAntiMatter == true)){

	fAntiProton_px = Track->Px();
	fAntiProton_py = Track->Py();
	fAntiProton_pz = Track->Pz();
	fAntiProton_px_Generated = Generated_px;
	fAntiProton_py_Generated = Generated_py;
	fAntiProton_pz_Generated = Generated_pz;
	fAntiProton_pTPC = Track->GetTPCmomentum();
	fAntiProton_Eta = Track->Eta();
	fAntiProton_Phi = Track->Phi();
	fAntiProton_TPC_Chi2 = Track->GetTPCchi2();
	fAntiProton_TPC_dEdx = Track->GetTPCsignal();
	if(fIsMC == false) fAntiProton_TPC_dEdx_Sigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
	if(fIsMC == true) fAntiProton_TPC_dEdx_Sigma = CalculateSigmadEdxTPC(*Track,3);
	fAntiProton_TOF_Mass2 = (Float_t)TOF_m2;
	fAntiProton_TOF_Mass2_Sigma = (Float_t)TOF_m2_Sigma;
	fAntiProton_ITS_dEdx = (Float_t)ITS_dEdx;
	fAntiProton_ITS_dEdx_Sigma = (Float_t)ITS_dEdx_Sigma;
	fAntiProton_DCAxy = DCAxy;
	fAntiProton_DCAz = DCAz;
	fAntiProton_TPC_nCrossedRows = Track->GetTPCCrossedRows();
	fAntiProton_TPC_nFindableCluster = Track->GetTPCNclsF();
	fAntiProton_TPC_nCluster = Track->GetTPCNcls();
	fAntiProton_ITS_nCluster = (UShort_t)ITS_nCluster;
	fAntiProton_PDG = PDG;
	fAntiProton_MotherPDG = MotherPDG;
	fAntiProton_FilterBit = Track->TestFilterBit(BIT(0));
	fAntiProton_ITS_Layer0 = Track->HasPointOnITSLayer(0);
	fAntiProton_ITS_Layer1 = Track->HasPointOnITSLayer(1);
	fAntiProton_ITS_Layer2 = Track->HasPointOnITSLayer(2);
	fAntiProton_ITS_Layer3 = Track->HasPointOnITSLayer(3);
	fAntiProton_ITS_Layer4 = Track->HasPointOnITSLayer(4);
	fAntiProton_ITS_Layer5 = Track->HasPointOnITSLayer(5);
	fAntiProton_Event_Multiplicity = Multiplicity;
	fAntiProton_Event_Centrality = Centrality;
	fAntiProton_Event_PrimaryVertexZ = PrimaryVertexZ;
	fAntiProton_Event_BField = IsPositiveBFieldPolarity;
	fAntiProton_Event_Identifier = EventID;
	fAntiProton_Event_RunNumber = RunNumber;
	fAntiProton_Event_IsFirstParticle = IsFirstAntiProton;
	fAntiProton_Event_TimeStamp = TimeStamp;
	fAntiProton_Event_RandomCrossCheckNumber = RandomCrossCheckNumber;
	fAntiProton_ID = track;

	fSaveTree_AntiProton->Fill();

	IsFirstAntiProton = false;

      } // end of isAntiProtons


    } // end of loop over Protons / AntiProtons



  } // end of SaveMatter / SaveAntiMatter







  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiProton);
  PostData(4,fSaveTree_AntiDeuteron);

} // end of UserExec















void AliAnalysisTask_pd_CreateTrees_PairsOnly::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
Double_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateBetaTOF(AliAODTrack &track)
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
Double_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateMassSquareTOF(AliAODTrack &track)
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








Double_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateSigmaMassSquareTOF(Double_t pT, Double_t massSq, Int_t ParticleSpecies, Int_t RunNumber)
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




















// apply track cuts for deuterons and antideuterons
Bool_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, Int_t ParticleSpecies = 0, Int_t RunNumber = 0)
{

  Bool_t PassedParticleCuts = false;


  // #######################################
  // ######## define particle cuts #########
  // #######################################

  Float_t Proton_pT_min = 0.0;
  Float_t Proton_pT_max = 3.0;
  Float_t Proton_eta_min = -0.8;
  Float_t Proton_eta_max = +0.8;
  Float_t Proton_DCAxy_max = 2.0; // cm
  Float_t Proton_DCAz_max = 0.1; // cm
  Float_t Proton_p_max = 30.0;
  Float_t Proton_TPC_RatioRowsFindableCluster_min = 0.8;
  Float_t Proton_TPC_dEdx_Sigma_max = 4.0;
  Float_t Proton_TPC_Chi2perCluster_max = 5.0;
  Float_t Proton_TPC_Chi2perNDF_max = 4.0;
  Int_t Proton_TPC_nCluster_min = 70;
  Int_t Proton_TPC_nCrossedRows_min = 60;
  Int_t Proton_TPC_nSharedCluster_max = 0;
  Float_t Proton_TPC_Threshold = 1.0;
  Float_t Proton_TOF_Mass2_Sigma_max = 4.0;
  Float_t Proton_ITS_dEdx_Sigma_max = 4.0;
  Int_t Proton_ITS_nCluster_min = 0;
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

  if(TMath::Abs(px) > Proton_p_max) return PassedParticleCuts;
  if(TMath::Abs(py) > Proton_p_max) return PassedParticleCuts;
  if(TMath::Abs(pz) > Proton_p_max) return PassedParticleCuts;



  // check if TOF information is available
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  // apply TOF-is-available cut (above threshold)
  if((pTPC >= Proton_TPC_Threshold) && (TOFisOK == false)) return PassedParticleCuts;


  // apply TPC Sigma cut
  if(fIsMC == false)  TPC_dEdx_Sigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
  if(fIsMC == true)   TPC_dEdx_Sigma = CalculateSigmadEdxTPC(Track,ParticleSpecies);
  if(TMath::IsNaN(TPC_dEdx_Sigma)) return PassedParticleCuts;

  if(fUseOpenCuts == false){
    if(TMath::Abs(TPC_dEdx_Sigma) > Proton_TPC_dEdx_Sigma_max) return PassedParticleCuts;
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
  if(TMath::Abs(DCAxy) > Proton_DCAxy_max) return PassedParticleCuts;

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Proton_DCAz_max) return PassedParticleCuts;

  // apply pT cut
  if(pT < Proton_pT_min || pT > Proton_pT_max) return PassedParticleCuts;

  Charge = Track.Charge();
  if(!(TMath::Abs(Charge) == 1)) return PassedParticleCuts;

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
  if((TMath::Abs(TPC_nCrossedRows) > 170) || (TPC_nCrossedRows < 0)) return PassedParticleCuts;
  if(TPC_nCrossedRows < Proton_TPC_nCrossedRows_min) return PassedParticleCuts;

  // apply zero shared cluster cut for TPC
  TPC_nSharedCluster = Track.GetTPCnclsS();
  if((TMath::Abs(TPC_nSharedCluster) > 10) || (TPC_nSharedCluster < 0)) return PassedParticleCuts;
  if(TPC_nSharedCluster > Proton_TPC_nSharedCluster_max) return PassedParticleCuts;

  // apply findable cluster cut for TPC
  TPC_nFindableCluster = Track.GetTPCNclsF();
  if((TMath::Abs(TPC_nFindableCluster) > 170) || (TPC_nFindableCluster < 0)) return PassedParticleCuts;
  if(TPC_nFindableCluster > 0) TPC_RatioRowsFindableCluster = ((Float_t)TPC_nCrossedRows / (Float_t)TPC_nFindableCluster);
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









// apply track cuts for deuterons and antideuterons
Bool_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, Int_t ParticleSpecies = 0, Int_t RunNumber = 0)
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
  if(fIsMC == false)  TPC_dEdx_Sigma = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
  if(fIsMC == true)   TPC_dEdx_Sigma = CalculateSigmadEdxTPC(Track,ParticleSpecies);

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








Double_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateSigmadEdxTPC(AliAODTrack &Track, Int_t ParticleSpecies){

  Double_t SigmaParticle = -999.0;
  Double_t SignalTPC = Track.GetTPCsignal();
  if(TMath::IsNaN(SignalTPC)) return SigmaParticle;

  Bool_t LHC20g7a   = false;
  Bool_t LHC20g7b   = false;
  Bool_t LHC22f3    = false;

  if(fIsMC == true){

    if(fCollisionSystem == 1) LHC20g7a	= true;
    if(fCollisionSystem == 2) LHC20g7b	= true;
    if(fCollisionSystem > 2)  LHC22f3	= true;

  } // end of fIsMC == true


  Bool_t isProton	= false;
  Bool_t isDeuteron     = false;
  Bool_t isAntiProton   = false;
  Bool_t isAntiDeuteron = false;

  if(ParticleSpecies == 1) isProton = true;
  if(ParticleSpecies == 2) isDeuteron = true;
  if(ParticleSpecies == 3) isAntiProton = true;
  if(ParticleSpecies == 4) isAntiDeuteron = true;


  Double_t pTPC = Track.GetTPCmomentum();
  if(TMath::IsNaN(pTPC)) return SigmaParticle;
  if(pTPC < 0.1) return SigmaParticle;
  if(pTPC > 6.0) return SigmaParticle;

  Double_t Mass = 0.0;
  if((isProton == true) || (isAntiProton == true))	Mass = AliPID::ParticleMass(AliPID::kProton);
  if((isDeuteron == true) || (isAntiDeuteron == true))	Mass = AliPID::ParticleMass(AliPID::kDeuteron);


  Double_t Resolution = 0.0;

  TF1 *Mean = new TF1("Mean","[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.0);
  Mean->FixParameter(5,1);
  Mean->FixParameter(6,Mass);


  // LHC20g7a
  if((LHC20g7a == true) && (isProton == true)){

    Mean->FixParameter(0,0.109964);
    Mean->FixParameter(1,362.867);
    Mean->FixParameter(2,5.49057e-14);
    Mean->FixParameter(3,2.06802);
    Mean->FixParameter(4,26.3114);
    Resolution = 0.0607035;


  }

  if((LHC20g7a == true) && (isDeuteron == true)){

    Mean->FixParameter(0,0.427825);
    Mean->FixParameter(1,73.1944);
    Mean->FixParameter(2,4.51003e-06);
    Mean->FixParameter(3,2.68815);
    Mean->FixParameter(4,28.7193);
    Resolution = 0.0584012;


  }

  if((LHC20g7a == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.158314);
    Mean->FixParameter(1,284.093);
    Mean->FixParameter(2,1.71436e-08);
    Mean->FixParameter(3,1.78061);
    Mean->FixParameter(4,-15.5739);
    Resolution = 0.059659;

  }

  if((LHC20g7a == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,0.474816);
    Mean->FixParameter(1,66.3673);
    Mean->FixParameter(2,3.26403e-06);
    Mean->FixParameter(3,2.70129);
    Mean->FixParameter(4,26.7397);
    Resolution = 0.0570914;


  }



  // LHC20g7b
  if((LHC20g7b == true) && (isProton == true)){

    Mean->FixParameter(0,0.243998);
    Mean->FixParameter(1,164.773);
    Mean->FixParameter(2,5.49057e-14);
    Mean->FixParameter(3,2.06802);
    Mean->FixParameter(4,26.3114);
    Resolution = 0.0629222;

  }

  if((LHC20g7b == true) && (isDeuteron == true)){

    Mean->FixParameter(0,0.58614);
    Mean->FixParameter(1,69.7817);
    Mean->FixParameter(2,1.99584);
    Mean->FixParameter(3,2.20193);
    Mean->FixParameter(4,15.2308);
    Resolution = 0.0544027;

  }

  if((LHC20g7b == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.241109);
    Mean->FixParameter(1,166.613);
    Mean->FixParameter(2,5.49057e-14);
    Mean->FixParameter(3,2.06802);
    Mean->FixParameter(4,26.3114);
    Resolution = 0.0551698;

  }

  if((LHC20g7b == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,0.428428);
    Mean->FixParameter(1,98.4929);
    Mean->FixParameter(2,4.18861);
    Mean->FixParameter(3,2.10266);
    Mean->FixParameter(4,15.5034);
    Resolution = 0.0551698;

  }






  // LHC22f3
  if((LHC22f3 == true) && (isProton == true)){

    Mean->FixParameter(0,0.279756);
    Mean->FixParameter(1,163.309);
    Mean->FixParameter(2,0.0273036);
    Mean->FixParameter(3,2.06802);
    Mean->FixParameter(4,26.3114);
    Resolution = 0.0705152;

  }


  if((LHC22f3 == true) && (isDeuteron == true)){

    Mean->FixParameter(0,0.474816);
    Mean->FixParameter(1,70.879);
    Mean->FixParameter(2,3.26403e-06);
    Mean->FixParameter(3,2.70129);
    Mean->FixParameter(4,26.7397);
    Resolution = 0.0651183;

  }

  if((LHC22f3 == true) && (isAntiProton == true)){

    Mean->FixParameter(0,0.187362);
    Mean->FixParameter(1,234.322);
    Mean->FixParameter(2,5.49057e-14);
    Mean->FixParameter(3,2.06802);
    Mean->FixParameter(4,26.3114);
    Resolution = 0.071852;

  }

  if((LHC22f3 == true) && (isAntiDeuteron == true)){

    Mean->FixParameter(0,0.474816);
    Mean->FixParameter(1,70.6882);
    Mean->FixParameter(2,3.26403e-06);
    Mean->FixParameter(3,2.70129);
    Mean->FixParameter(4,26.7397);
    Resolution = 0.0646782;

  }



  Double_t mean = Mean->Eval(pTPC);
  Mean->Delete();



  Double_t ScaleFactor = 1.0-(Resolution);
  Double_t sigma = (mean*ScaleFactor) - mean;
  if(TMath::Abs(sigma) < 0.0001) return -999.0;

  SigmaParticle = (mean - SignalTPC) / (sigma);
  if(TMath::IsNaN(SigmaParticle)) return -999.0;

  return SigmaParticle;


} // end of CalculateSigmadEdxTPC





Double_t AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateSigmadEdxITS(AliAODTrack &Track, Int_t ParticleSpecies, Int_t RunNumber){

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








