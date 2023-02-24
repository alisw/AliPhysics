#include "TChain.h" 
#include "TTree.h"
#include "TList.h"
#include "TH2F.h"
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
  fProton_TPC_Chi2(0),
  fProton_TPC_dEdx(0),
  fProton_TPC_dEdx_nSigma(0),
  fProton_TOF_Beta(0),
  fProton_TOF_Beta_nSigma(0),
  fProton_TOF_Mass2(0),
  fProton_TOF_Mass2_nSigma(0),
  fProton_ITS_dEdx(0),
  fProton_ITS_dEdx_nSigma(0),
  fProton_DCAxy(0),
  fProton_DCAz(0),
  fProton_Event_Centrality(0),
  fProton_Event_PrimaryVertexZ(0),
  fProton_Event_BField(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nSharedCluster(0),
  fProton_TPC_nFindableCluster(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_Event_nParticles(0),
  fProton_ID(0),
  fProton_Event_Identifier(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(0),
  fDeuteron_py(0),
  fDeuteron_pz(0),
  fDeuteron_pTPC(0),
  fDeuteron_Eta(0),
  fDeuteron_Phi(0),
  fDeuteron_TPC_Chi2(0),
  fDeuteron_TPC_dEdx(0),
  fDeuteron_TPC_dEdx_nSigma(0),
  fDeuteron_TOF_Beta(0),
  fDeuteron_TOF_Beta_nSigma(0),
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
  fDeuteron_Event_nParticles(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Identifier(0),
  fSaveTree_AntiProton(0),
  fAntiProton_px(0),
  fAntiProton_py(0),
  fAntiProton_pz(0),
  fAntiProton_pTPC(0),
  fAntiProton_Eta(0),
  fAntiProton_Phi(0),
  fAntiProton_TPC_Chi2(0),
  fAntiProton_TPC_dEdx(0),
  fAntiProton_TPC_dEdx_nSigma(0),
  fAntiProton_TOF_Beta(0),
  fAntiProton_TOF_Beta_nSigma(0),
  fAntiProton_TOF_Mass2(0),
  fAntiProton_TOF_Mass2_nSigma(0),
  fAntiProton_ITS_dEdx(0),
  fAntiProton_ITS_dEdx_nSigma(0),
  fAntiProton_DCAxy(0),
  fAntiProton_DCAz(0),
  fAntiProton_Event_Centrality(0),
  fAntiProton_Event_PrimaryVertexZ(0),
  fAntiProton_Event_BField(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nSharedCluster(0),
  fAntiProton_TPC_nFindableCluster(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_Event_nParticles(0),
  fAntiProton_ID(0),
  fAntiProton_Event_Identifier(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(0),
  fAntiDeuteron_py(0),
  fAntiDeuteron_pz(0),
  fAntiDeuteron_pTPC(0),
  fAntiDeuteron_Eta(0),
  fAntiDeuteron_Phi(0),
  fAntiDeuteron_TPC_Chi2(0),
  fAntiDeuteron_TPC_dEdx(0),
  fAntiDeuteron_TPC_dEdx_nSigma(0),
  fAntiDeuteron_TOF_Beta(0),
  fAntiDeuteron_TOF_Beta_nSigma(0),
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
  fAntiDeuteron_Event_nParticles(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Identifier(0),
  fHistoList(0),
  h_Proton_TOF_m2_NoTOFcut(0),
  h_Deuteron_TOF_m2_NoTOFcut(0),
  h_AntiProton_TOF_m2_NoTOFcut(0),
  h_AntiDeuteron_TOF_m2_NoTOFcut(0),
  h_Proton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut(0)
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
  fProton_TPC_Chi2(0),
  fProton_TPC_dEdx(0),
  fProton_TPC_dEdx_nSigma(0),
  fProton_TOF_Beta(0),
  fProton_TOF_Beta_nSigma(0),
  fProton_TOF_Mass2(0),
  fProton_TOF_Mass2_nSigma(0),
  fProton_ITS_dEdx(0),
  fProton_ITS_dEdx_nSigma(0),
  fProton_DCAxy(0),
  fProton_DCAz(0),
  fProton_Event_Centrality(0),
  fProton_Event_PrimaryVertexZ(0),
  fProton_Event_BField(0),
  fProton_TPC_nCrossedRows(0),
  fProton_TPC_nSharedCluster(0),
  fProton_TPC_nFindableCluster(0),
  fProton_TPC_nCluster(0),
  fProton_ITS_nCluster(0),
  fProton_Event_nParticles(0),
  fProton_ID(0),
  fProton_Event_Identifier(0),
  fSaveTree_Deuteron(0),
  fDeuteron_px(0),
  fDeuteron_py(0),
  fDeuteron_pz(0),
  fDeuteron_pTPC(0),
  fDeuteron_Eta(0),
  fDeuteron_Phi(0),
  fDeuteron_TPC_Chi2(0),
  fDeuteron_TPC_dEdx(0),
  fDeuteron_TPC_dEdx_nSigma(0),
  fDeuteron_TOF_Beta(0),
  fDeuteron_TOF_Beta_nSigma(0),
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
  fDeuteron_Event_nParticles(0),
  fDeuteron_ID(0),
  fDeuteron_Event_Identifier(0),
  fSaveTree_AntiProton(0),
  fAntiProton_px(0),
  fAntiProton_py(0),
  fAntiProton_pz(0),
  fAntiProton_pTPC(0),
  fAntiProton_Eta(0),
  fAntiProton_Phi(0),
  fAntiProton_TPC_Chi2(0),
  fAntiProton_TPC_dEdx(0),
  fAntiProton_TPC_dEdx_nSigma(0),
  fAntiProton_TOF_Beta(0),
  fAntiProton_TOF_Beta_nSigma(0),
  fAntiProton_TOF_Mass2(0),
  fAntiProton_TOF_Mass2_nSigma(0),
  fAntiProton_ITS_dEdx(0),
  fAntiProton_ITS_dEdx_nSigma(0),
  fAntiProton_DCAxy(0),
  fAntiProton_DCAz(0),
  fAntiProton_Event_Centrality(0),
  fAntiProton_Event_PrimaryVertexZ(0),
  fAntiProton_Event_BField(0),
  fAntiProton_TPC_nCrossedRows(0),
  fAntiProton_TPC_nSharedCluster(0),
  fAntiProton_TPC_nFindableCluster(0),
  fAntiProton_TPC_nCluster(0),
  fAntiProton_ITS_nCluster(0),
  fAntiProton_Event_nParticles(0),
  fAntiProton_ID(0),
  fAntiProton_Event_Identifier(0),
  fSaveTree_AntiDeuteron(0),
  fAntiDeuteron_px(0),
  fAntiDeuteron_py(0),
  fAntiDeuteron_pz(0),
  fAntiDeuteron_pTPC(0),
  fAntiDeuteron_Eta(0),
  fAntiDeuteron_Phi(0),
  fAntiDeuteron_TPC_Chi2(0),
  fAntiDeuteron_TPC_dEdx(0),
  fAntiDeuteron_TPC_dEdx_nSigma(0),
  fAntiDeuteron_TOF_Beta(0),
  fAntiDeuteron_TOF_Beta_nSigma(0),
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
  fAntiDeuteron_Event_nParticles(0),
  fAntiDeuteron_ID(0),
  fAntiDeuteron_Event_Identifier(0),
  fHistoList(0),
  h_Proton_TOF_m2_NoTOFcut(0),
  h_Deuteron_TOF_m2_NoTOFcut(0),
  h_AntiProton_TOF_m2_NoTOFcut(0),
  h_AntiDeuteron_TOF_m2_NoTOFcut(0),
  h_Proton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_Deuteron_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());
  DefineOutput(5,TList::Class());

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

  if(fHistoList)
    {
      delete fHistoList;
    }

}







void AliAnalysisTask_pd_CreateTrees_PairsOnly::UserCreateOutputObjects()
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


  fSaveTree_Proton = new TTree("fSaveTree_Proton","fSaveTree_Proton");
  fSaveTree_Proton->Branch("Proton_px",&fProton_px,"Proton_px/f");
  fSaveTree_Proton->Branch("Proton_py",&fProton_py,"Proton_py/f");
  fSaveTree_Proton->Branch("Proton_pz",&fProton_pz,"Proton_pz/f");
  fSaveTree_Proton->Branch("Proton_pTPC",&fProton_pTPC,"Proton_pTPC/f");
  fSaveTree_Proton->Branch("Proton_Eta",&fProton_Eta,"Proton_Eta/f");
  fSaveTree_Proton->Branch("Proton_Phi",&fProton_Phi,"Proton_Phi/f");
  fSaveTree_Proton->Branch("Proton_TPC_Chi2",&fProton_TPC_Chi2,"Proton_TPC_Chi2/f");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx",&fProton_TPC_dEdx,"Proton_TPC_dEdx/f");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx_nSigma",&fProton_TPC_dEdx_nSigma,"Proton_TPC_dEdx_nSigma/f");
  fSaveTree_Proton->Branch("Proton_TOF_Beta",&fProton_TOF_Beta,"Proton_TOF_Beta/f");
  fSaveTree_Proton->Branch("Proton_TOF_Beta_nSigma",&fProton_TOF_Beta_nSigma,"Proton_TOF_Beta_nSigma/f");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2",&fProton_TOF_Mass2,"Proton_TOF_Mass2/f");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2_nSigma",&fProton_TOF_Mass2_nSigma,"Proton_TOF_Mass2_nSigma/f");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx",&fProton_ITS_dEdx,"Proton_ITS_dEdx/f");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx_nSigma",&fProton_ITS_dEdx_nSigma,"Proton_ITS_dEdx_nSigma/f");
  fSaveTree_Proton->Branch("Proton_DCAxy",&fProton_DCAxy,"Proton_DCAxy/f");
  fSaveTree_Proton->Branch("Proton_DCAz",&fProton_DCAz,"Proton_DCAz/f");
  fSaveTree_Proton->Branch("Proton_Event_Centrality",&fProton_Event_Centrality,"Proton_Event_Centrality/f");
  fSaveTree_Proton->Branch("Proton_Event_PrimaryVertexZ",&fProton_Event_PrimaryVertexZ,"Proton_Event_PrimaryVertexZ/f");
  fSaveTree_Proton->Branch("Proton_Event_BField",&fProton_Event_BField,"Proton_Event_BField/f");
  fSaveTree_Proton->Branch("Proton_TPC_nCrossedRows",&fProton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fSaveTree_Proton->Branch("Proton_TPC_nSharedCluster",&fProton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fSaveTree_Proton->Branch("Proton_TPC_nFindableCluster",&fProton_TPC_nFindableCluster,"Proton_TPC_nFindableCluster/s");
  fSaveTree_Proton->Branch("Proton_TPC_nCluster",&fProton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fSaveTree_Proton->Branch("Proton_ITS_nCluster",&fProton_ITS_nCluster,"Proton_ITS_nCluster/S");
  fSaveTree_Proton->Branch("Proton_Event_nParticles",&fProton_Event_nParticles,"Proton_Event_nParticles/s");
  fSaveTree_Proton->Branch("Proton_ID",&fProton_ID,"Proton_ID/i");
  fSaveTree_Proton->Branch("Proton_Event_Identifier",&fProton_Event_Identifier,"Proton_Event_Identifier/i");


  fSaveTree_Deuteron = new TTree("fSaveTree_Deuteron","fSaveTree_Deuteron");
  fSaveTree_Deuteron->Branch("Deuteron_px",&fDeuteron_px,"Deuteron_px/f");
  fSaveTree_Deuteron->Branch("Deuteron_py",&fDeuteron_py,"Deuteron_py/f");
  fSaveTree_Deuteron->Branch("Deuteron_pz",&fDeuteron_pz,"Deuteron_pz/f");
  fSaveTree_Deuteron->Branch("Deuteron_pTPC",&fDeuteron_pTPC,"Deuteron_pTPC/f");
  fSaveTree_Deuteron->Branch("Deuteron_Eta",&fDeuteron_Eta,"Deuteron_Eta/f");
  fSaveTree_Deuteron->Branch("Deuteron_Phi",&fDeuteron_Phi,"Deuteron_Phi/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_Chi2",&fDeuteron_TPC_Chi2,"Deuteron_TPC_Chi2/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx",&fDeuteron_TPC_dEdx,"Deuteron_TPC_dEdx/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_dEdx_nSigma",&fDeuteron_TPC_dEdx_nSigma,"Deuteron_TPC_dEdx_nSigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Beta",&fDeuteron_TOF_Beta,"Deuteron_TOF_Beta/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Beta_nSigma",&fDeuteron_TOF_Beta_nSigma,"Deuteron_TOF_Beta_nSigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2",&fDeuteron_TOF_Mass2,"Deuteron_TOF_Mass2/f");
  fSaveTree_Deuteron->Branch("Deuteron_TOF_Mass2_nSigma",&fDeuteron_TOF_Mass2_nSigma,"Deuteron_TOF_Mass2_nSigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx",&fDeuteron_ITS_dEdx,"Deuteron_ITS_dEdx/f");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_dEdx_nSigma",&fDeuteron_ITS_dEdx_nSigma,"Deuteron_ITS_dEdx_nSigma/f");
  fSaveTree_Deuteron->Branch("Deuteron_DCAxy",&fDeuteron_DCAxy,"Deuteron_DCAxy/f");
  fSaveTree_Deuteron->Branch("Deuteron_DCAz",&fDeuteron_DCAz,"Deuteron_DCAz/f");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Centrality",&fDeuteron_Event_Centrality,"Deuteron_Event_Centrality/f");
  fSaveTree_Deuteron->Branch("Deuteron_Event_PrimaryVertexZ",&fDeuteron_Event_PrimaryVertexZ,"Deuteron_Event_PrimaryVertexZ/f");
  fSaveTree_Deuteron->Branch("Deuteron_Event_BField",&fDeuteron_Event_BField,"Deuteron_Event_BField/f");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&fDeuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&fDeuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nFindableCluster",&fDeuteron_TPC_nFindableCluster,"Deuteron_TPC_nFindableCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_TPC_nCluster",&fDeuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fSaveTree_Deuteron->Branch("Deuteron_ITS_nCluster",&fDeuteron_ITS_nCluster,"Deuteron_ITS_nCluster/S");
  fSaveTree_Deuteron->Branch("Deuteron_Event_nParticles",&fDeuteron_Event_nParticles,"Deuteron_Event_nParticles/s");
  fSaveTree_Deuteron->Branch("Deuteron_ID",&fDeuteron_ID,"Deuteron_ID/i");
  fSaveTree_Deuteron->Branch("Deuteron_Event_Identifier",&fDeuteron_Event_Identifier,"Deuteron_Event_Identifier/i");




  fSaveTree_AntiProton = new TTree("fSaveTree_AntiProton","fSaveTree_AntiProton");
  fSaveTree_AntiProton->Branch("AntiProton_px",&fAntiProton_px,"AntiProton_px/f");
  fSaveTree_AntiProton->Branch("AntiProton_py",&fAntiProton_py,"AntiProton_py/f");
  fSaveTree_AntiProton->Branch("AntiProton_pz",&fAntiProton_pz,"AntiProton_pz/f");
  fSaveTree_AntiProton->Branch("AntiProton_pTPC",&fAntiProton_pTPC,"AntiProton_pTPC/f");
  fSaveTree_AntiProton->Branch("AntiProton_Eta",&fAntiProton_Eta,"AntiProton_Eta/f");
  fSaveTree_AntiProton->Branch("AntiProton_Phi",&fAntiProton_Phi,"AntiProton_Phi/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_Chi2",&fAntiProton_TPC_Chi2,"AntiProton_TPC_Chi2/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx",&fAntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx_nSigma",&fAntiProton_TPC_dEdx_nSigma,"AntiProton_TPC_dEdx_nSigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Beta",&fAntiProton_TOF_Beta,"AntiProton_TOF_Beta/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Beta_nSigma",&fAntiProton_TOF_Beta_nSigma,"AntiProton_TOF_Beta_nSigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2",&fAntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/f");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2_nSigma",&fAntiProton_TOF_Mass2_nSigma,"AntiProton_TOF_Mass2_nSigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx",&fAntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/f");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx_nSigma",&fAntiProton_ITS_dEdx_nSigma,"AntiProton_ITS_dEdx_nSigma/f");
  fSaveTree_AntiProton->Branch("AntiProton_DCAxy",&fAntiProton_DCAxy,"AntiProton_DCAxy/f");
  fSaveTree_AntiProton->Branch("AntiProton_DCAz",&fAntiProton_DCAz,"AntiProton_DCAz/f");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Centrality",&fAntiProton_Event_Centrality,"AntiProton_Event_Centrality/f");
  fSaveTree_AntiProton->Branch("AntiProton_Event_PrimaryVertexZ",&fAntiProton_Event_PrimaryVertexZ,"AntiProton_Event_PrimaryVertexZ/f");
  fSaveTree_AntiProton->Branch("AntiProton_Event_BField",&fAntiProton_Event_BField,"AntiProton_Event_BField/f");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&fAntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&fAntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nFindableCluster",&fAntiProton_TPC_nFindableCluster,"AntiProton_TPC_nFindableCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCluster",&fAntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_nCluster",&fAntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/S");
  fSaveTree_AntiProton->Branch("AntiProton_Event_nParticles",&fAntiProton_Event_nParticles,"AntiProton_Event_nParticles/s");
  fSaveTree_AntiProton->Branch("AntiProton_ID",&fAntiProton_ID,"AntiProton_ID/i");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Identifier",&fAntiProton_Event_Identifier,"AntiProton_Event_Identifier/i");


  fSaveTree_AntiDeuteron = new TTree("fSaveTree_AntiDeuteron","fSaveTree_AntiDeuteron");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_px",&fAntiDeuteron_px,"AntiDeuteron_px/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_py",&fAntiDeuteron_py,"AntiDeuteron_py/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pz",&fAntiDeuteron_pz,"AntiDeuteron_pz/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&fAntiDeuteron_pTPC,"AntiDeuteron_pTPC/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&fAntiDeuteron_Eta,"AntiDeuteron_Eta/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&fAntiDeuteron_Phi,"AntiDeuteron_Phi/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2",&fAntiDeuteron_TPC_Chi2,"AntiDeuteron_TPC_Chi2/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&fAntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_nSigma",&fAntiDeuteron_TPC_dEdx_nSigma,"AntiDeuteron_TPC_dEdx_nSigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta",&fAntiDeuteron_TOF_Beta,"AntiDeuteron_TOF_Beta/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta_nSigma",&fAntiDeuteron_TOF_Beta_nSigma,"AntiDeuteron_TOF_Beta_nSigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&fAntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_nSigma",&fAntiDeuteron_TOF_Mass2_nSigma,"AntiDeuteron_TOF_Mass2_nSigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&fAntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_nSigma",&fAntiDeuteron_ITS_dEdx_nSigma,"AntiDeuteron_ITS_dEdx_nSigma/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&fAntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&fAntiDeuteron_DCAz,"AntiDeuteron_DCAz/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Centrality",&fAntiDeuteron_Event_Centrality,"AntiDeuteron_Event_Centrality/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_PrimaryVertexZ",&fAntiDeuteron_Event_PrimaryVertexZ,"AntiDeuteron_Event_PrimaryVertexZ/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_BField",&fAntiDeuteron_Event_BField,"AntiDeuteron_Event_BField/f");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&fAntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&fAntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nFindableCluster",&fAntiDeuteron_TPC_nFindableCluster,"AntiDeuteron_TPC_nFindableCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&fAntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&fAntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/S");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_nParticles",&fAntiDeuteron_Event_nParticles,"AntiDeuteron_Event_nParticles/s");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_ID",&fAntiDeuteron_ID,"AntiDeuteron_ID/i");
  fSaveTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&fAntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/i");




  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiProton);
  PostData(4,fSaveTree_AntiDeuteron);
  PostData(5,fHistoList);



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


  // debug analysis
  bool DebugEventSelection  = false;

 

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



  // get primary vertex
  AliAODVertex *PrimaryVertex = fAODEvent->GetPrimaryVertex();
  if(!PrimaryVertex)::Warning("AliAnalsisTask_pd_CreateTrees_PairsOnlyd::UserExec","No AliAODVertex object found!");
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
  if((Centrality < Centrality_min) || (Centrality > Centrality_max)) return;




  // get event information
  unsigned int PeriodNumber	= fAODEvent->GetPeriodNumber();
  unsigned int OrbitNumber	= fAODEvent->GetOrbitNumber();
  unsigned int BunchCrossNumber	= fAODEvent->GetBunchCrossNumber();
  int RunNumber			= fAODEvent->GetRunNumber();
  float BField			= fAODEvent->GetMagneticField();

  if(TMath::IsNaN(BField)) return;



  // print event information
  if(DebugEventSelection)
  {

    cout << "" << endl;
    cout << "fCollisionSystem:\t\t" << fCollisionSystem << std::endl;
    cout << "PeriodNumber:\t\t\t" << PeriodNumber << endl;
    cout << "OrbitNumber:\t\t\t" << OrbitNumber << endl;
    cout << "BunchCrossNumber:\t\t" << BunchCrossNumber << endl;
    cout << "Centrality / Multiplicity:\t" << Centrality << " %" << endl;
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
  float     Proton_TPC_Chi2;
  float     Proton_TPC_dEdx;
  float     Proton_TPC_dEdx_nSigma;
  float     Proton_TOF_Beta;
  float     Proton_TOF_Beta_nSigma;
  float     Proton_TOF_Mass2;
  float     Proton_TOF_Mass2_nSigma;
  float     Proton_ITS_dEdx;
  float     Proton_ITS_dEdx_nSigma;
  float     Proton_DCAxy;
  float     Proton_DCAz;
  unsigned short    Proton_TPC_nCrossedRows;
  unsigned short    Proton_TPC_nSharedCluster;
  unsigned short    Proton_TPC_nFindableCluster;
  unsigned short    Proton_TPC_nCluster;
  short    Proton_ITS_nCluster;
  unsigned int      Proton_ID;
  unsigned int      Proton_Event_Identifier;


  TTree *fTempTree_Proton = new TTree("fTempTree_Proton","fTempTree_Proton");
  fTempTree_Proton->Branch("Proton_px",&Proton_px,"Proton_px/f");
  fTempTree_Proton->Branch("Proton_py",&Proton_py,"Proton_py/f");
  fTempTree_Proton->Branch("Proton_pz",&Proton_pz,"Proton_pz/f");
  fTempTree_Proton->Branch("Proton_pTPC",&Proton_pTPC,"Proton_pTPC/f");
  fTempTree_Proton->Branch("Proton_Eta",&Proton_Eta,"Proton_Eta/f");
  fTempTree_Proton->Branch("Proton_Phi",&Proton_Phi,"Proton_Phi/f");
  fTempTree_Proton->Branch("Proton_TPC_Chi2",&Proton_TPC_Chi2,"Proton_TPC_Chi2/f");
  fTempTree_Proton->Branch("Proton_TPC_dEdx",&Proton_TPC_dEdx,"Proton_TPC_dEdx/f");
  fTempTree_Proton->Branch("Proton_TPC_dEdx_nSigma",&Proton_TPC_dEdx_nSigma,"Proton_TPC_dEdx_nSigma/f");
  fTempTree_Proton->Branch("Proton_TOF_Beta",&Proton_TOF_Beta,"Proton_TOF_Beta/f");
  fTempTree_Proton->Branch("Proton_TOF_Beta_nSigma",&Proton_TOF_Beta_nSigma,"Proton_TOF_Beta_nSigma/f");
  fTempTree_Proton->Branch("Proton_TOF_Mass2",&Proton_TOF_Mass2,"Proton_TOF_Mass2/f");
  fTempTree_Proton->Branch("Proton_TOF_Mass2_nSigma",&Proton_TOF_Mass2_nSigma,"Proton_TOF_Mass2_nSigma/f");
  fTempTree_Proton->Branch("Proton_ITS_dEdx",&Proton_ITS_dEdx,"Proton_ITS_dEdx/f");
  fTempTree_Proton->Branch("Proton_ITS_dEdx_nSigma",&Proton_ITS_dEdx_nSigma,"Proton_ITS_dEdx_nSigma/f");
  fTempTree_Proton->Branch("Proton_DCAxy",&Proton_DCAxy,"Proton_DCAxy/f");
  fTempTree_Proton->Branch("Proton_DCAz",&Proton_DCAz,"Proton_DCAz/f");
  fTempTree_Proton->Branch("Proton_TPC_nCrossedRows",&Proton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fTempTree_Proton->Branch("Proton_TPC_nSharedCluster",&Proton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fTempTree_Proton->Branch("Proton_TPC_nFindableCluster",&Proton_TPC_nFindableCluster,"Proton_TPC_nFindableCluster/s");
  fTempTree_Proton->Branch("Proton_TPC_nCluster",&Proton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fTempTree_Proton->Branch("Proton_ITS_nCluster",&Proton_ITS_nCluster,"Proton_ITS_nCluster/S");
  fTempTree_Proton->Branch("Proton_ID",&Proton_ID,"Proton_ID/i");
  fTempTree_Proton->Branch("Proton_Event_Identifier",&Proton_Event_Identifier,"Proton_Event_Identifier/i");

  unsigned short nProtonsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// proton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;

    // apply proton cuts
    bool PassedProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,true,RunNumber);
    if(!PassedProtonCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    double TOF_Beta	    = -999.0;
    double TOF_Beta_nSigma  = -999.0;
    double TOF_m2	    = -999.0;
    double TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    bool TOFisOK = false;
    if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
    if(TOFisOK){

      TOF_Beta	      = CalculateBetaTOF(*Track);
      TOF_Beta_nSigma = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
      TOF_m2	      = CalculateMassSquareTOF(*Track);
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,1,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = -999;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,1,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();

    }


    Proton_px			    = Track->Px();
    Proton_py			    = Track->Py();
    Proton_pz			    = Track->Pz();
    Proton_pTPC			    = Track->GetTPCmomentum();
    Proton_Eta			    = Track->Eta();
    Proton_Phi			    = Track->Phi();
    Proton_TPC_Chi2		    = Track->GetTPCchi2();
    Proton_TPC_dEdx		    = Track->GetTPCsignal();
    Proton_TPC_dEdx_nSigma	    = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    Proton_TOF_Beta		    = TOF_Beta;
    Proton_TOF_Beta_nSigma	    = TOF_Beta_nSigma;
    Proton_TOF_Mass2		    = TOF_m2;
    Proton_TOF_Mass2_nSigma	    = TOF_m2_nSigma;
    Proton_ITS_dEdx		    = ITS_dEdx;
    Proton_ITS_dEdx_nSigma	    = ITS_dEdx_nSigma;
    Proton_DCAxy		    = DCAxy;
    Proton_DCAz			    = DCAz;
    Proton_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    Proton_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    Proton_TPC_nFindableCluster	    = Track->GetTPCNclsF();
    Proton_TPC_nCluster		    = Track->GetTPCNcls();
    Proton_ITS_nCluster		    = ITS_nCluster;
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
  float     Deuteron_TPC_Chi2;
  float     Deuteron_TPC_dEdx;
  float     Deuteron_TPC_dEdx_nSigma;
  float     Deuteron_TOF_Beta;
  float     Deuteron_TOF_Beta_nSigma;
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
  short  Deuteron_ITS_nCluster;
  unsigned int    Deuteron_ID;
  unsigned int    Deuteron_Event_Identifier;

  TTree *fTempTree_Deuteron = new TTree("fTempTree_Deuteron","fTempTree_Deuteron");
  fTempTree_Deuteron->Branch("Deuteron_px",&Deuteron_px,"Deuteron_px/f");
  fTempTree_Deuteron->Branch("Deuteron_py",&Deuteron_py,"Deuteron_py/f");
  fTempTree_Deuteron->Branch("Deuteron_pz",&Deuteron_pz,"Deuteron_pz/f");
  fTempTree_Deuteron->Branch("Deuteron_pTPC",&Deuteron_pTPC,"Deuteron_pTPC/f");
  fTempTree_Deuteron->Branch("Deuteron_Eta",&Deuteron_Eta,"Deuteron_Eta/f");
  fTempTree_Deuteron->Branch("Deuteron_Phi",&Deuteron_Phi,"Deuteron_Phi/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_Chi2",&Deuteron_TPC_Chi2,"Deuteron_TPC_Chi2/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_dEdx",&Deuteron_TPC_dEdx,"Deuteron_TPC_dEdx/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_dEdx_nSigma",&Deuteron_TPC_dEdx_nSigma,"Deuteron_TPC_dEdx_nSigma/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Beta",&Deuteron_TOF_Beta,"Deuteron_TOF_Beta/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Beta_nSigma",&Deuteron_TOF_Beta_nSigma,"Deuteron_TOF_Beta_nSigma/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Mass2",&Deuteron_TOF_Mass2,"Deuteron_TOF_Mass2/f");
  fTempTree_Deuteron->Branch("Deuteron_TOF_Mass2_nSigma",&Deuteron_TOF_Mass2_nSigma,"Deuteron_TOF_Mass2_nSigma/f");
  fTempTree_Deuteron->Branch("Deuteron_ITS_dEdx",&Deuteron_ITS_dEdx,"Deuteron_ITS_dEdx/f");
  fTempTree_Deuteron->Branch("Deuteron_ITS_dEdx_nSigma",&Deuteron_ITS_dEdx_nSigma,"Deuteron_ITS_dEdx_nSigma/f");
  fTempTree_Deuteron->Branch("Deuteron_DCAxy",&Deuteron_DCAxy,"Deuteron_DCAxy/f");
  fTempTree_Deuteron->Branch("Deuteron_DCAz",&Deuteron_DCAz,"Deuteron_DCAz/f");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nCrossedRows",&Deuteron_TPC_nCrossedRows,"Deuteron_TPC_nCrossedRows/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nSharedCluster",&Deuteron_TPC_nSharedCluster,"Deuteron_TPC_nSharedCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nFindableCluster",&Deuteron_TPC_nFindableCluster,"Deuteron_TPC_nFindableCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_TPC_nCluster",&Deuteron_TPC_nCluster,"Deuteron_TPC_nCluster/s");
  fTempTree_Deuteron->Branch("Deuteron_ITS_nCluster",&Deuteron_ITS_nCluster,"Deuteron_ITS_nCluster/S");
  fTempTree_Deuteron->Branch("Deuteron_ID",&Deuteron_ID,"Deuteron_ID/i");
  fTempTree_Deuteron->Branch("Deuteron_Event_Identifier",&Deuteron_Event_Identifier,"Deuteron_Event_Identifier/i");

  unsigned short nDeuteronsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// deuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply deuteron cuts
    bool PassedDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,true,RunNumber);
    if(!PassedDeuteronCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    double TOF_Beta	    = -999.0;
    double TOF_Beta_nSigma  = -999.0;
    double TOF_m2	    = -999.0;
    double TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    bool TOFisOK = false;
    if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
    if(TOFisOK){

      TOF_Beta	      = CalculateBetaTOF(*Track);
      TOF_Beta_nSigma = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
      TOF_m2	      = CalculateMassSquareTOF(*Track);
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,2,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = -999;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,2,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();

    }


    Deuteron_px			    = Track->Px();
    Deuteron_py			    = Track->Py();
    Deuteron_pz			    = Track->Pz();
    Deuteron_pTPC		    = Track->GetTPCmomentum();
    Deuteron_Eta		    = Track->Eta();
    Deuteron_Phi		    = Track->Phi();
    Deuteron_TPC_Chi2		    = Track->GetTPCchi2();
    Deuteron_TPC_dEdx		    = Track->GetTPCsignal();
    Deuteron_TPC_dEdx_nSigma	    = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    Deuteron_TOF_Beta		    = TOF_Beta;
    Deuteron_TOF_Beta_nSigma	    = TOF_Beta_nSigma;
    Deuteron_TOF_Mass2		    = TOF_m2;
    Deuteron_TOF_Mass2_nSigma	    = TOF_m2_nSigma;
    Deuteron_ITS_dEdx		    = ITS_dEdx;
    Deuteron_ITS_dEdx_nSigma	    = ITS_dEdx_nSigma;
    Deuteron_DCAxy		    = DCAxy;
    Deuteron_DCAz		    = DCAz;
    Deuteron_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    Deuteron_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    Deuteron_TPC_nFindableCluster   = Track->GetTPCNclsF();
    Deuteron_TPC_nCluster	    = Track->GetTPCNcls();
    Deuteron_ITS_nCluster	    = ITS_nCluster;
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
      TBranch *Branch_Proton_TPC_Chi2		  = fTempTree_Proton->GetBranch("Proton_TPC_Chi2");
      TBranch *Branch_Proton_TPC_dEdx		  = fTempTree_Proton->GetBranch("Proton_TPC_dEdx");
      TBranch *Branch_Proton_TPC_dEdx_nSigma	  = fTempTree_Proton->GetBranch("Proton_TPC_dEdx_nSigma");
      TBranch *Branch_Proton_TOF_Beta		  = fTempTree_Proton->GetBranch("Proton_TOF_Beta");
      TBranch *Branch_Proton_TOF_Beta_nSigma	  = fTempTree_Proton->GetBranch("Proton_TOF_Beta_nSigma");
      TBranch *Branch_Proton_TOF_Mass2		  = fTempTree_Proton->GetBranch("Proton_TOF_Mass2");
      TBranch *Branch_Proton_TOF_Mass2_nSigma	  = fTempTree_Proton->GetBranch("Proton_TOF_Mass2_nSigma");
      TBranch *Branch_Proton_ITS_dEdx		  = fTempTree_Proton->GetBranch("Proton_ITS_dEdx");
      TBranch *Branch_Proton_ITS_dEdx_nSigma	  = fTempTree_Proton->GetBranch("Proton_ITS_dEdx_nSigma");
      TBranch *Branch_Proton_DCAxy		  = fTempTree_Proton->GetBranch("Proton_DCAxy");
      TBranch *Branch_Proton_DCAz		  = fTempTree_Proton->GetBranch("Proton_DCAz");
      TBranch *Branch_Proton_TPC_nCrossedRows	  = fTempTree_Proton->GetBranch("Proton_TPC_nCrossedRows");
      TBranch *Branch_Proton_TPC_nSharedCluster	  = fTempTree_Proton->GetBranch("Proton_TPC_nSharedCluster");
      TBranch *Branch_Proton_TPC_nFindableCluster = fTempTree_Proton->GetBranch("Proton_TPC_nFindableCluster");
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
      Branch_Proton_TPC_Chi2->SetAddress(&fProton_TPC_Chi2);
      Branch_Proton_TPC_dEdx->SetAddress(&fProton_TPC_dEdx);
      Branch_Proton_TPC_dEdx_nSigma->SetAddress(&fProton_TPC_dEdx_nSigma);
      Branch_Proton_TOF_Beta->SetAddress(&fProton_TOF_Beta);
      Branch_Proton_TOF_Beta_nSigma->SetAddress(&fProton_TOF_Beta_nSigma);
      Branch_Proton_TOF_Mass2->SetAddress(&fProton_TOF_Mass2);
      Branch_Proton_TOF_Mass2_nSigma->SetAddress(&fProton_TOF_Mass2_nSigma);
      Branch_Proton_ITS_dEdx->SetAddress(&fProton_ITS_dEdx);
      Branch_Proton_ITS_dEdx_nSigma->SetAddress(&fProton_ITS_dEdx_nSigma);
      Branch_Proton_DCAxy->SetAddress(&fProton_DCAxy);
      Branch_Proton_DCAz->SetAddress(&fProton_DCAz);
      Branch_Proton_TPC_nCrossedRows->SetAddress(&fProton_TPC_nCrossedRows);
      Branch_Proton_TPC_nSharedCluster->SetAddress(&fProton_TPC_nSharedCluster);
      Branch_Proton_TPC_nFindableCluster->SetAddress(&fProton_TPC_nFindableCluster);
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
      Branch_Proton_TPC_Chi2->SetAutoDelete(true);
      Branch_Proton_TPC_dEdx->SetAutoDelete(true);
      Branch_Proton_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_Proton_TOF_Beta->SetAutoDelete(true);
      Branch_Proton_TOF_Beta_nSigma->SetAutoDelete(true);
      Branch_Proton_TOF_Mass2->SetAutoDelete(true);
      Branch_Proton_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_Proton_ITS_dEdx->SetAutoDelete(true);
      Branch_Proton_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_Proton_DCAxy->SetAutoDelete(true);
      Branch_Proton_DCAz->SetAutoDelete(true);
      Branch_Proton_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_Proton_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_Proton_TPC_nFindableCluster->SetAutoDelete(true);
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
      Branch_Proton_TPC_Chi2->GetEntry(Proton);
      Branch_Proton_TPC_dEdx->GetEntry(Proton);
      Branch_Proton_TPC_dEdx_nSigma->GetEntry(Proton);
      Branch_Proton_TOF_Beta->GetEntry(Proton);
      Branch_Proton_TOF_Beta_nSigma->GetEntry(Proton);
      Branch_Proton_TOF_Mass2->GetEntry(Proton);
      Branch_Proton_TOF_Mass2_nSigma->GetEntry(Proton);
      Branch_Proton_ITS_dEdx->GetEntry(Proton);
      Branch_Proton_ITS_dEdx_nSigma->GetEntry(Proton);
      Branch_Proton_DCAxy->GetEntry(Proton);
      Branch_Proton_DCAz->GetEntry(Proton);
      Branch_Proton_TPC_nCrossedRows->GetEntry(Proton);
      Branch_Proton_TPC_nSharedCluster->GetEntry(Proton);
      Branch_Proton_TPC_nFindableCluster->GetEntry(Proton);
      Branch_Proton_TPC_nCluster->GetEntry(Proton);
      Branch_Proton_ITS_nCluster->GetEntry(Proton);
      Branch_Proton_ID->GetEntry(Proton);
      Branch_Proton_Event_Identifier->GetEntry(Proton);

      fProton_Event_nParticles	    = nProtonsSelected;
      fProton_Event_Centrality	    = Centrality;
      fProton_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fProton_Event_BField	    = BField;
      

      fSaveTree_Proton->Fill();

    } // end of loop (copy protons)


    for(int Deuteron = 0; Deuteron < nDeuteronsSelected; Deuteron++){

      TBranch *Branch_Deuteron_px		    = fTempTree_Deuteron->GetBranch("Deuteron_px");
      TBranch *Branch_Deuteron_py		    = fTempTree_Deuteron->GetBranch("Deuteron_py");
      TBranch *Branch_Deuteron_pz		    = fTempTree_Deuteron->GetBranch("Deuteron_pz");
      TBranch *Branch_Deuteron_pTPC		    = fTempTree_Deuteron->GetBranch("Deuteron_pTPC");
      TBranch *Branch_Deuteron_Eta		    = fTempTree_Deuteron->GetBranch("Deuteron_Eta");
      TBranch *Branch_Deuteron_Phi		    = fTempTree_Deuteron->GetBranch("Deuteron_Phi");
      TBranch *Branch_Deuteron_TPC_Chi2		    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_Chi2");
      TBranch *Branch_Deuteron_TPC_dEdx		    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_dEdx");
      TBranch *Branch_Deuteron_TPC_dEdx_nSigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TPC_dEdx_nSigma");
      TBranch *Branch_Deuteron_TOF_Beta		    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Beta");
      TBranch *Branch_Deuteron_TOF_Beta_nSigma	    = fTempTree_Deuteron->GetBranch("Deuteron_TOF_Beta_nSigma");
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
      TBranch *Branch_Deuteron_ID		    = fTempTree_Deuteron->GetBranch("Deuteron_ID");
      TBranch *Branch_Deuteron_Event_Identifier	    = fTempTree_Deuteron->GetBranch("Deuteron_Event_Identifier");


      Branch_Deuteron_px->SetAddress(&fDeuteron_px);
      Branch_Deuteron_py->SetAddress(&fDeuteron_py);
      Branch_Deuteron_pz->SetAddress(&fDeuteron_pz);
      Branch_Deuteron_pTPC->SetAddress(&fDeuteron_pTPC);
      Branch_Deuteron_Eta->SetAddress(&fDeuteron_Eta);
      Branch_Deuteron_Phi->SetAddress(&fDeuteron_Phi);
      Branch_Deuteron_TPC_Chi2->SetAddress(&fDeuteron_TPC_Chi2);
      Branch_Deuteron_TPC_dEdx->SetAddress(&fDeuteron_TPC_dEdx);
      Branch_Deuteron_TPC_dEdx_nSigma->SetAddress(&fDeuteron_TPC_dEdx_nSigma);
      Branch_Deuteron_TOF_Beta->SetAddress(&fDeuteron_TOF_Beta);
      Branch_Deuteron_TOF_Beta_nSigma->SetAddress(&fDeuteron_TOF_Beta_nSigma);
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
      Branch_Deuteron_TOF_Beta->SetAutoDelete(true);
      Branch_Deuteron_TOF_Beta_nSigma->SetAutoDelete(true);
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
      Branch_Deuteron_ID->SetAutoDelete(true);
      Branch_Deuteron_Event_Identifier->SetAutoDelete(true);

      Branch_Deuteron_px->GetEntry(Deuteron);
      Branch_Deuteron_py->GetEntry(Deuteron);
      Branch_Deuteron_pz->GetEntry(Deuteron);
      Branch_Deuteron_pTPC->GetEntry(Deuteron);
      Branch_Deuteron_Eta->GetEntry(Deuteron);
      Branch_Deuteron_Phi->GetEntry(Deuteron);
      Branch_Deuteron_TPC_Chi2->GetEntry(Deuteron);
      Branch_Deuteron_TPC_dEdx->GetEntry(Deuteron);
      Branch_Deuteron_TPC_dEdx_nSigma->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Beta->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Beta_nSigma->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Mass2->GetEntry(Deuteron);
      Branch_Deuteron_TOF_Mass2_nSigma->GetEntry(Deuteron);
      Branch_Deuteron_ITS_dEdx->GetEntry(Deuteron);
      Branch_Deuteron_ITS_dEdx_nSigma->GetEntry(Deuteron);
      Branch_Deuteron_DCAxy->GetEntry(Deuteron);
      Branch_Deuteron_DCAz->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nCrossedRows->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nSharedCluster->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nFindableCluster->GetEntry(Deuteron);
      Branch_Deuteron_TPC_nCluster->GetEntry(Deuteron);
      Branch_Deuteron_ITS_nCluster->GetEntry(Deuteron);
      Branch_Deuteron_ID->GetEntry(Deuteron);
      Branch_Deuteron_Event_Identifier->GetEntry(Deuteron);

      fDeuteron_Event_nParticles      = nDeuteronsSelected;
      fDeuteron_Event_Centrality      = Centrality;
      fDeuteron_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fDeuteron_Event_BField	      = BField;

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
  float     AntiProton_TPC_Chi2;
  float     AntiProton_TPC_dEdx;
  float     AntiProton_TPC_dEdx_nSigma;
  float     AntiProton_TOF_Beta;
  float     AntiProton_TOF_Beta_nSigma;
  float     AntiProton_TOF_Mass2;
  float     AntiProton_TOF_Mass2_nSigma;
  float     AntiProton_ITS_dEdx;
  float     AntiProton_ITS_dEdx_nSigma;
  float     AntiProton_DCAxy;
  float     AntiProton_DCAz;
  unsigned short    AntiProton_TPC_nCrossedRows;
  unsigned short    AntiProton_TPC_nSharedCluster;
  unsigned short    AntiProton_TPC_nFindableCluster;
  unsigned short    AntiProton_TPC_nCluster;
  short    AntiProton_ITS_nCluster;
  unsigned int      AntiProton_ID;
  unsigned int      AntiProton_Event_Identifier;

  fTempTree_AntiProton = new TTree("fTempTree_AntiProton","fTempTree_AntiProton");
  fTempTree_AntiProton->Branch("AntiProton_px",&AntiProton_px,"AntiProton_px/f");
  fTempTree_AntiProton->Branch("AntiProton_py",&AntiProton_py,"AntiProton_py/f");
  fTempTree_AntiProton->Branch("AntiProton_pz",&AntiProton_pz,"AntiProton_pz/f");
  fTempTree_AntiProton->Branch("AntiProton_pTPC",&AntiProton_pTPC,"AntiProton_pTPC/f");
  fTempTree_AntiProton->Branch("AntiProton_Eta",&AntiProton_Eta,"AntiProton_Eta/f");
  fTempTree_AntiProton->Branch("AntiProton_Phi",&AntiProton_Phi,"AntiProton_Phi/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_Chi2",&AntiProton_TPC_Chi2,"AntiProton_TPC_Chi2/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_dEdx",&AntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_dEdx_nSigma",&AntiProton_TPC_dEdx_nSigma,"AntiProton_TPC_dEdx_nSigma/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Beta",&AntiProton_TOF_Beta,"AntiProton_TOF_Beta/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Beta_nSigma",&AntiProton_TOF_Beta_nSigma,"AntiProton_TOF_Beta_nSigma/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Mass2",&AntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/f");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Mass2_nSigma",&AntiProton_TOF_Mass2_nSigma,"AntiProton_TOF_Mass2_nSigma/f");
  fTempTree_AntiProton->Branch("AntiProton_ITS_dEdx",&AntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/f");
  fTempTree_AntiProton->Branch("AntiProton_ITS_dEdx_nSigma",&AntiProton_ITS_dEdx_nSigma,"AntiProton_ITS_dEdx_nSigma/f");
  fTempTree_AntiProton->Branch("AntiProton_DCAxy",&AntiProton_DCAxy,"AntiProton_DCAxy/f");
  fTempTree_AntiProton->Branch("AntiProton_DCAz",&AntiProton_DCAz,"AntiProton_DCAz/f");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&AntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&AntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nFindableCluster",&AntiProton_TPC_nFindableCluster,"AntiProton_TPC_nFindableCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nCluster",&AntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_ITS_nCluster",&AntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/S");
  fTempTree_AntiProton->Branch("AntiProton_ID",&AntiProton_ID,"AntiProton_ID/i");
  fTempTree_AntiProton->Branch("AntiProton_Event_Identifier",&AntiProton_Event_Identifier,"AntiProton_Event_Identifier/i");

  unsigned short nAntiProtonsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// antiproton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antiproton cuts
    bool PassedAntiProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,false,RunNumber);
    if(!PassedAntiProtonCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    double TOF_Beta	    = -999.0;
    double TOF_Beta_nSigma  = -999.0;
    double TOF_m2	    = -999.0;
    double TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    bool TOFisOK = false;
    if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
    if(TOFisOK){

      TOF_Beta	      = CalculateBetaTOF(*Track);
      TOF_Beta_nSigma = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
      TOF_m2	      = CalculateMassSquareTOF(*Track);
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,3,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = -999;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,3,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();

    }


    AntiProton_px		      = Track->Px();
    AntiProton_py		      = Track->Py();
    AntiProton_pz		      = Track->Pz();
    AntiProton_pTPC		      = Track->GetTPCmomentum();
    AntiProton_Eta		      = Track->Eta();
    AntiProton_Phi		      = Track->Phi();
    AntiProton_TPC_Chi2		      = Track->GetTPCchi2();
    AntiProton_TPC_dEdx		      = Track->GetTPCsignal();
    AntiProton_TPC_dEdx_nSigma	      = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    AntiProton_TOF_Beta		      = TOF_Beta;
    AntiProton_TOF_Beta_nSigma	      = TOF_Beta_nSigma;
    AntiProton_TOF_Mass2	      = TOF_m2;
    AntiProton_TOF_Mass2_nSigma	      = TOF_m2_nSigma;
    AntiProton_ITS_dEdx		      = ITS_dEdx;
    AntiProton_ITS_dEdx_nSigma	      = ITS_dEdx_nSigma;
    AntiProton_DCAxy		      = DCAxy;
    AntiProton_DCAz		      = DCAz;
    AntiProton_TPC_nCrossedRows	      = Track->GetTPCCrossedRows();
    AntiProton_TPC_nSharedCluster     = Track->GetTPCnclsS();
    AntiProton_TPC_nFindableCluster   = Track->GetTPCNclsF();
    AntiProton_TPC_nCluster	      = Track->GetTPCNcls();
    AntiProton_ITS_nCluster	      = ITS_nCluster;
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
    float     AntiDeuteron_TPC_Chi2;
    float     AntiDeuteron_TPC_dEdx;
    float     AntiDeuteron_TPC_dEdx_nSigma;
    float     AntiDeuteron_TOF_Beta;
    float     AntiDeuteron_TOF_Beta_nSigma;
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
    short  AntiDeuteron_ITS_nCluster;
    unsigned int    AntiDeuteron_ID;
    unsigned int    AntiDeuteron_Event_Identifier;

  fTempTree_AntiDeuteron = new TTree("fTempTree_AntiDeuteron","fTempTree_AntiDeuteron");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_px",&AntiDeuteron_px,"AntiDeuteron_px/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_py",&AntiDeuteron_py,"AntiDeuteron_py/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pz",&AntiDeuteron_pz,"AntiDeuteron_pz/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_pTPC",&AntiDeuteron_pTPC,"AntiDeuteron_pTPC/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Eta",&AntiDeuteron_Eta,"AntiDeuteron_Eta/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Phi",&AntiDeuteron_Phi,"AntiDeuteron_Phi/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_Chi2",&AntiDeuteron_TPC_Chi2,"AntiDeuteron_TPC_Chi2/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx",&AntiDeuteron_TPC_dEdx,"AntiDeuteron_TPC_dEdx/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_dEdx_nSigma",&AntiDeuteron_TPC_dEdx_nSigma,"AntiDeuteron_TPC_dEdx_nSigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta",&AntiDeuteron_TOF_Beta,"AntiDeuteron_TOF_Beta/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Beta_nSigma",&AntiDeuteron_TOF_Beta_nSigma,"AntiDeuteron_TOF_Beta_nSigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2",&AntiDeuteron_TOF_Mass2,"AntiDeuteron_TOF_Mass2/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TOF_Mass2_nSigma",&AntiDeuteron_TOF_Mass2_nSigma,"AntiDeuteron_TOF_Mass2_nSigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx",&AntiDeuteron_ITS_dEdx,"AntiDeuteron_ITS_dEdx/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_dEdx_nSigma",&AntiDeuteron_ITS_dEdx_nSigma,"AntiDeuteron_ITS_dEdx_nSigma/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_DCAxy",&AntiDeuteron_DCAxy,"AntiDeuteron_DCAxy/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_DCAz",&AntiDeuteron_DCAz,"AntiDeuteron_DCAz/f");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCrossedRows",&AntiDeuteron_TPC_nCrossedRows,"AntiDeuteron_TPC_nCrossedRows/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nSharedCluster",&AntiDeuteron_TPC_nSharedCluster,"AntiDeuteron_TPC_nSharedCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nFindableCluster",&AntiDeuteron_TPC_nFindableCluster,"AntiDeuteron_TPC_nFindableCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_TPC_nCluster",&AntiDeuteron_TPC_nCluster,"AntiDeuteron_TPC_nCluster/s");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ITS_nCluster",&AntiDeuteron_ITS_nCluster,"AntiDeuteron_ITS_nCluster/S");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_ID",&AntiDeuteron_ID,"AntiDeuteron_ID/i");
  fTempTree_AntiDeuteron->Branch("AntiDeuteron_Event_Identifier",&AntiDeuteron_Event_Identifier,"AntiDeuteron_Event_Identifier/i");




  unsigned short nAntiDeuteronsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// antideuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pd_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antideuteron cuts
    bool PassedAntiDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,false,RunNumber);
    if(!PassedAntiDeuteronCuts) continue;
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    double TOF_Beta	    = -999.0;
    double TOF_Beta_nSigma  = -999.0;
    double TOF_m2	    = -999.0;
    double TOF_m2_nSigma    = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    bool TOFisOK = false;
    if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;
    
    if(TOFisOK){

      TOF_Beta	      = CalculateBetaTOF(*Track);
      TOF_Beta_nSigma = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
      TOF_m2	      = CalculateMassSquareTOF(*Track);
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,4,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = -999;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,4,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();

    }


    AntiDeuteron_px			= Track->Px();
    AntiDeuteron_py			= Track->Py();
    AntiDeuteron_pz			= Track->Pz();
    AntiDeuteron_pTPC			= Track->GetTPCmomentum();
    AntiDeuteron_Eta			= Track->Eta();
    AntiDeuteron_Phi			= Track->Phi();
    AntiDeuteron_TPC_Chi2		= Track->GetTPCchi2();
    AntiDeuteron_TPC_dEdx		= Track->GetTPCsignal();
    AntiDeuteron_TPC_dEdx_nSigma	= fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron);
    AntiDeuteron_TOF_Beta		= TOF_Beta;
    AntiDeuteron_TOF_Beta_nSigma	= TOF_Beta_nSigma;
    AntiDeuteron_TOF_Mass2		= TOF_m2;
    AntiDeuteron_TOF_Mass2_nSigma	= TOF_m2_nSigma;
    AntiDeuteron_ITS_dEdx		= ITS_dEdx;
    AntiDeuteron_ITS_dEdx_nSigma	= ITS_dEdx_nSigma;
    AntiDeuteron_DCAxy			= DCAxy;
    AntiDeuteron_DCAz			= DCAz;
    AntiDeuteron_TPC_nCrossedRows	= Track->GetTPCCrossedRows();
    AntiDeuteron_TPC_nSharedCluster	= Track->GetTPCnclsS();
    AntiDeuteron_TPC_nFindableCluster	= Track->GetTPCNclsF();
    AntiDeuteron_TPC_nCluster		= Track->GetTPCNcls();
    AntiDeuteron_ITS_nCluster		= ITS_nCluster;
    AntiDeuteron_ID			= Track->GetID();
    AntiDeuteron_Event_Identifier	= BunchCrossNumber + (OrbitNumber*3564) + (PeriodNumber*16777215*3564);
 
    fTempTree_AntiDeuteron->Fill();
    nAntiDeuteronsSelected++;

  } // end of antideuteron loop




  if((nAntiProtonsSelected > 0) && (nAntiDeuteronsSelected > 0)){

    for(int AntiProton = 0; AntiProton < nAntiProtonsSelected; AntiProton++){

      TBranch *Branch_AntiProton_px		      = fTempTree_AntiProton->GetBranch("AntiProton_px");
      TBranch *Branch_AntiProton_py		      = fTempTree_AntiProton->GetBranch("AntiProton_py");
      TBranch *Branch_AntiProton_pz		      = fTempTree_AntiProton->GetBranch("AntiProton_pz");
      TBranch *Branch_AntiProton_pTPC		      = fTempTree_AntiProton->GetBranch("AntiProton_pTPC");
      TBranch *Branch_AntiProton_Eta		      = fTempTree_AntiProton->GetBranch("AntiProton_Eta");
      TBranch *Branch_AntiProton_Phi		      = fTempTree_AntiProton->GetBranch("AntiProton_Phi");
      TBranch *Branch_AntiProton_TPC_Chi2	      = fTempTree_AntiProton->GetBranch("AntiProton_TPC_Chi2");
      TBranch *Branch_AntiProton_TPC_dEdx	      = fTempTree_AntiProton->GetBranch("AntiProton_TPC_dEdx");
      TBranch *Branch_AntiProton_TPC_dEdx_nSigma      = fTempTree_AntiProton->GetBranch("AntiProton_TPC_dEdx_nSigma");
      TBranch *Branch_AntiProton_TOF_Beta	      = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Beta");
      TBranch *Branch_AntiProton_TOF_Beta_nSigma      = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Beta_nSigma");
      TBranch *Branch_AntiProton_TOF_Mass2	      = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Mass2");
      TBranch *Branch_AntiProton_TOF_Mass2_nSigma     = fTempTree_AntiProton->GetBranch("AntiProton_TOF_Mass2_nSigma");
      TBranch *Branch_AntiProton_ITS_dEdx	      = fTempTree_AntiProton->GetBranch("AntiProton_ITS_dEdx");
      TBranch *Branch_AntiProton_ITS_dEdx_nSigma      = fTempTree_AntiProton->GetBranch("AntiProton_ITS_dEdx_nSigma");
      TBranch *Branch_AntiProton_DCAxy		      = fTempTree_AntiProton->GetBranch("AntiProton_DCAxy");
      TBranch *Branch_AntiProton_DCAz		      = fTempTree_AntiProton->GetBranch("AntiProton_DCAz");
      TBranch *Branch_AntiProton_TPC_nCrossedRows     = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nCrossedRows");
      TBranch *Branch_AntiProton_TPC_nSharedCluster   = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nSharedCluster");
      TBranch *Branch_AntiProton_TPC_nFindableCluster = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nFindableCluster");
      TBranch *Branch_AntiProton_TPC_nCluster	      = fTempTree_AntiProton->GetBranch("AntiProton_TPC_nCluster");
      TBranch *Branch_AntiProton_ITS_nCluster	      = fTempTree_AntiProton->GetBranch("AntiProton_ITS_nCluster");
      TBranch *Branch_AntiProton_ID		      = fTempTree_AntiProton->GetBranch("AntiProton_ID");
      TBranch *Branch_AntiProton_Event_Identifier     = fTempTree_AntiProton->GetBranch("AntiProton_Event_Identifier");


      Branch_AntiProton_px->SetAddress(&fAntiProton_px);
      Branch_AntiProton_py->SetAddress(&fAntiProton_py);
      Branch_AntiProton_pz->SetAddress(&fAntiProton_pz);
      Branch_AntiProton_pTPC->SetAddress(&fAntiProton_pTPC);
      Branch_AntiProton_Eta->SetAddress(&fAntiProton_Eta);
      Branch_AntiProton_Phi->SetAddress(&fAntiProton_Phi);
      Branch_AntiProton_TPC_Chi2->SetAddress(&fAntiProton_TPC_Chi2);
      Branch_AntiProton_TPC_dEdx->SetAddress(&fAntiProton_TPC_dEdx);
      Branch_AntiProton_TPC_dEdx_nSigma->SetAddress(&fAntiProton_TPC_dEdx_nSigma);
      Branch_AntiProton_TOF_Beta->SetAddress(&fAntiProton_TOF_Beta);
      Branch_AntiProton_TOF_Beta_nSigma->SetAddress(&fAntiProton_TOF_Beta_nSigma);
      Branch_AntiProton_TOF_Mass2->SetAddress(&fAntiProton_TOF_Mass2);
      Branch_AntiProton_TOF_Mass2_nSigma->SetAddress(&fAntiProton_TOF_Mass2_nSigma);
      Branch_AntiProton_ITS_dEdx->SetAddress(&fAntiProton_ITS_dEdx);
      Branch_AntiProton_ITS_dEdx_nSigma->SetAddress(&fAntiProton_ITS_dEdx_nSigma);
      Branch_AntiProton_DCAxy->SetAddress(&fAntiProton_DCAxy);
      Branch_AntiProton_DCAz->SetAddress(&fAntiProton_DCAz);
      Branch_AntiProton_TPC_nCrossedRows->SetAddress(&fAntiProton_TPC_nCrossedRows);
      Branch_AntiProton_TPC_nSharedCluster->SetAddress(&fAntiProton_TPC_nSharedCluster);
      Branch_AntiProton_TPC_nFindableCluster->SetAddress(&fAntiProton_TPC_nFindableCluster);
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
      Branch_AntiProton_TPC_Chi2->SetAutoDelete(true);
      Branch_AntiProton_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiProton_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiProton_TOF_Beta->SetAutoDelete(true);
      Branch_AntiProton_TOF_Beta_nSigma->SetAutoDelete(true);
      Branch_AntiProton_TOF_Mass2->SetAutoDelete(true);
      Branch_AntiProton_TOF_Mass2_nSigma->SetAutoDelete(true);
      Branch_AntiProton_ITS_dEdx->SetAutoDelete(true);
      Branch_AntiProton_ITS_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiProton_DCAxy->SetAutoDelete(true);
      Branch_AntiProton_DCAz->SetAutoDelete(true);
      Branch_AntiProton_TPC_nCrossedRows->SetAutoDelete(true);
      Branch_AntiProton_TPC_nSharedCluster->SetAutoDelete(true);
      Branch_AntiProton_TPC_nFindableCluster->SetAutoDelete(true);
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
      Branch_AntiProton_TPC_Chi2->GetEntry(AntiProton);
      Branch_AntiProton_TPC_dEdx->GetEntry(AntiProton);
      Branch_AntiProton_TPC_dEdx_nSigma->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Beta->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Beta_nSigma->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Mass2->GetEntry(AntiProton);
      Branch_AntiProton_TOF_Mass2_nSigma->GetEntry(AntiProton);
      Branch_AntiProton_ITS_dEdx->GetEntry(AntiProton);
      Branch_AntiProton_ITS_dEdx_nSigma->GetEntry(AntiProton);
      Branch_AntiProton_DCAxy->GetEntry(AntiProton);
      Branch_AntiProton_DCAz->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nCrossedRows->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nSharedCluster->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nFindableCluster->GetEntry(AntiProton);
      Branch_AntiProton_TPC_nCluster->GetEntry(AntiProton);
      Branch_AntiProton_ITS_nCluster->GetEntry(AntiProton);
      Branch_AntiProton_ID->GetEntry(AntiProton);
      Branch_AntiProton_Event_Identifier->GetEntry(AntiProton);

      fAntiProton_Event_nParticles	= nAntiProtonsSelected;
      fAntiProton_Event_Centrality	= Centrality;
      fAntiProton_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fAntiProton_Event_BField		= BField;

      fSaveTree_AntiProton->Fill();

    } // end of loop (copy antiprotons)


    for(int AntiDeuteron = 0; AntiDeuteron < nAntiDeuteronsSelected; AntiDeuteron++){

      TBranch *Branch_AntiDeuteron_px			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_px");
      TBranch *Branch_AntiDeuteron_py			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_py");
      TBranch *Branch_AntiDeuteron_pz			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pz");
      TBranch *Branch_AntiDeuteron_pTPC			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_pTPC");
      TBranch *Branch_AntiDeuteron_Eta			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Eta");
      TBranch *Branch_AntiDeuteron_Phi			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Phi");
      TBranch *Branch_AntiDeuteron_TPC_Chi2		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_Chi2");
      TBranch *Branch_AntiDeuteron_TPC_dEdx		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_dEdx");
      TBranch *Branch_AntiDeuteron_TPC_dEdx_nSigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TPC_dEdx_nSigma");
      TBranch *Branch_AntiDeuteron_TOF_Beta		= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Beta");
      TBranch *Branch_AntiDeuteron_TOF_Beta_nSigma	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_TOF_Beta_nSigma");
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
      TBranch *Branch_AntiDeuteron_ID			= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_ID");
      TBranch *Branch_AntiDeuteron_Event_Identifier	= fTempTree_AntiDeuteron->GetBranch("AntiDeuteron_Event_Identifier");


      Branch_AntiDeuteron_px->SetAddress(&fAntiDeuteron_px);
      Branch_AntiDeuteron_py->SetAddress(&fAntiDeuteron_py);
      Branch_AntiDeuteron_pz->SetAddress(&fAntiDeuteron_pz);
      Branch_AntiDeuteron_pTPC->SetAddress(&fAntiDeuteron_pTPC);
      Branch_AntiDeuteron_Eta->SetAddress(&fAntiDeuteron_Eta);
      Branch_AntiDeuteron_Phi->SetAddress(&fAntiDeuteron_Phi);
      Branch_AntiDeuteron_TPC_Chi2->SetAddress(&fAntiDeuteron_TPC_Chi2);
      Branch_AntiDeuteron_TPC_dEdx->SetAddress(&fAntiDeuteron_TPC_dEdx);
      Branch_AntiDeuteron_TPC_dEdx_nSigma->SetAddress(&fAntiDeuteron_TPC_dEdx_nSigma);
      Branch_AntiDeuteron_TOF_Beta->SetAddress(&fAntiDeuteron_TOF_Beta);
      Branch_AntiDeuteron_TOF_Beta_nSigma->SetAddress(&fAntiDeuteron_TOF_Beta_nSigma);
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
      Branch_AntiDeuteron_ID->SetAddress(&fAntiDeuteron_ID);
      Branch_AntiDeuteron_Event_Identifier->SetAddress(&fAntiDeuteron_Event_Identifier);


      Branch_AntiDeuteron_px->SetAutoDelete(true);
      Branch_AntiDeuteron_py->SetAutoDelete(true);
      Branch_AntiDeuteron_pz->SetAutoDelete(true);
      Branch_AntiDeuteron_pTPC->SetAutoDelete(true);
      Branch_AntiDeuteron_Eta->SetAutoDelete(true);
      Branch_AntiDeuteron_Phi->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_Chi2->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_dEdx->SetAutoDelete(true);
      Branch_AntiDeuteron_TPC_dEdx_nSigma->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Beta->SetAutoDelete(true);
      Branch_AntiDeuteron_TOF_Beta_nSigma->SetAutoDelete(true);
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
      Branch_AntiDeuteron_ID->SetAutoDelete(true);
      Branch_AntiDeuteron_Event_Identifier->SetAutoDelete(true);

      Branch_AntiDeuteron_px->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_py->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_pz->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_pTPC->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Eta->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Phi->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_Chi2->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_dEdx->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_dEdx_nSigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Beta->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Beta_nSigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Mass2->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TOF_Mass2_nSigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ITS_dEdx->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ITS_dEdx_nSigma->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_DCAxy->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_DCAz->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nCrossedRows->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nSharedCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nFindableCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_TPC_nCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ITS_nCluster->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_ID->GetEntry(AntiDeuteron);
      Branch_AntiDeuteron_Event_Identifier->GetEntry(AntiDeuteron);

      fAntiDeuteron_Event_nParticles	  = nAntiDeuteronsSelected;
      fAntiDeuteron_Event_Centrality	  = Centrality;
      fAntiDeuteron_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fAntiDeuteron_Event_BField	  = BField;

      fSaveTree_AntiDeuteron->Fill();


    } // end of loop (copy antideuterons)


  } // end of same-event

  fTempTree_AntiProton->Delete();
  fTempTree_AntiDeuteron->Delete();












  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_Deuteron);
  PostData(3,fSaveTree_AntiProton);
  PostData(4,fSaveTree_AntiDeuteron);
  PostData(5,fHistoList);

} // end of UserExec















void AliAnalysisTask_pd_CreateTrees_PairsOnly::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateBetaTOF(AliAODTrack &track)
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
double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateMassSquareTOF(AliAODTrack &track)
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








double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateSigmaMassSquareTOF(double pT, double massSq, int ParticleSpecies, int RunNumber)
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

  if(ParticleSpecies == 1) isProton	  = true;
  if(ParticleSpecies == 2) isDeuteron	  = true;
  if(ParticleSpecies == 3) isAntiProton	  = true;
  if(ParticleSpecies == 4) isAntiDeuteron = true;
  

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

  // MetaLHC18q PbPb data 
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




  double mean = Mean->Eval(pT);
  double sigma = Sigma->Eval(pT);

  Mean->Delete();
  Sigma->Delete();


  SigmaParticle = (massSq - mean)/(sigma);
  return SigmaParticle;

} // end of CalculateSigmaMassSquareTOF










// apply track cuts for protons and antiprotons
bool AliAnalysisTask_pd_CreateTrees_PairsOnly::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber)
{

  bool PassedParticleCuts = false;

  // define deuteron and antideuteron track cuts
  double Proton_pT_min = 0.0;
  double Proton_pT_max = 4.0;
  double Proton_eta_min = -0.8;
  double Proton_eta_max = +0.8;
  double Proton_DCAxy_max = 0.2; // cm
  double Proton_DCAz_max = 0.1; // cm

  double Proton_TPC_RatioRowsFindableCluster_min = 0.83;
  double Proton_TPC_dEdx_nSigma_max = 3.0;
  double Proton_TPC_Chi2perCluster_max = 4.0;
  double Proton_TPC_Chi2perNDF_max = 4.0;
  int Proton_TPC_nCluster_min = 80;
  int Proton_TPC_nCrossedRows_min = 70;
  int Proton_TPC_nSharedCluster_max = 0;
  double Proton_TPC_Threshold = 0.7;

  double Proton_TOF_m2_nSigma_max = 3.0;
  double Proton_TOF_m2_nSigma_max_low_pTPC = 7.0;

  double Proton_ITS_dEdx_nSigma_max = 3.0;
  int Proton_ITS_nCluster_min = 2;

  bool UseTOF = true;
  bool UseITS = true;


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
  double DCAxy = xv[0];
  double DCAz = xv[1];
  if(TMath::IsNaN(DCAxy)) return PassedParticleCuts;
  if(TMath::IsNaN(DCAz)) return PassedParticleCuts;
  
  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Proton_DCAxy_max) return PassedParticleCuts;

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Proton_DCAz_max) return PassedParticleCuts;

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

































// apply track cuts for deuterons and antideuterons
bool AliAnalysisTask_pd_CreateTrees_PairsOnly::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber)
{

  bool PassedParticleCuts = false;

  // define deuteron and antideuteron track cuts
  double Deuteron_pT_min = 0.0;
  double Deuteron_pT_max = 4.0;
  double Deuteron_eta_min = -0.8;
  double Deuteron_eta_max = +0.8;
  double Deuteron_DCAxy_max = 0.2; // cm
  double Deuteron_DCAz_max = 0.1; // cm

  double Deuteron_TPC_RatioRowsFindableCluster_min = 0.83;
  double Deuteron_TPC_dEdx_nSigma_max = 3.0;
  double Deuteron_TPC_Chi2perCluster_max = 4.0;
  double Deuteron_TPC_Chi2perNDF_max = 4.0;
  int Deuteron_TPC_nCluster_min = 80;
  int Deuteron_TPC_nCrossedRows_min = 70;
  int Deuteron_TPC_nSharedCluster_max = 0;
  double Deuteron_TPC_Threshold = 1.0;

  double Deuteron_TOF_m2_nSigma_max = 3.0;
  double Deuteron_TOF_m2_nSigma_max_low_pTPC = 7.0;

  double Deuteron_ITS_dEdx_nSigma_max = 3.0;
  int Deuteron_ITS_nCluster_min = 2;

  bool UseTOF = true;
  bool UseITS = true;


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
  if(TMath::Abs(TPC_dEdx_nSigma) > Deuteron_TPC_dEdx_nSigma_max) return PassedParticleCuts;

  // get DCA information
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  double DCAxy = xv[0];
  double DCAz = xv[1];
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


    // apply tight TOF m2 cut above pTPC threshold   
    if(pTPC >= Deuteron_TPC_Threshold){

      if(TMath::Abs(TOF_m2_nSigma) > Deuteron_TOF_m2_nSigma_max) return PassedParticleCuts;

    }

    // apply loose TOF m2 cut below pTPC threshold
    if(pTPC < Deuteron_TPC_Threshold){

      if(TMath::Abs(TOF_m2_nSigma) > Deuteron_TOF_m2_nSigma_max_low_pTPC) return PassedParticleCuts;

    }
 
  } // end of TOFisOK













  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckDeuteronCuts











double AliAnalysisTask_pd_CreateTrees_PairsOnly::CalculateSigmadEdxITS(AliAODTrack &Track, int ParticleSpecies, int RunNumber){

  double SigmaParticle = -999.0;
  double SignalITS = Track.GetITSsignal();
  if(TMath::IsNaN(SignalITS)) return SigmaParticle;

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

  if(ParticleSpecies == 1) isProton = true;
  if(ParticleSpecies == 2) isDeuteron = true;
  if(ParticleSpecies == 3) isAntiProton = true;
  if(ParticleSpecies == 4) isAntiDeuteron = true;

  double p = Track.P();

  TF1 *Mean = new TF1("Mean","[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.01,6.0);

  if((isProton == true) || (isAntiProton == true)){

    Mean->FixParameter(0,2.36861e-07);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,9.55834);
    Mean->FixParameter(4,17081);
    Mean->FixParameter(5,1);
    Mean->FixParameter(6,0.93827208816);

  }

  if((isDeuteron == true) || (isAntiDeuteron == true)){

    Mean->FixParameter(0,7.41722e-06);
    Mean->FixParameter(1,-55831.1);
    Mean->FixParameter(2,-238672);
    Mean->FixParameter(3,11249.3);
    Mean->FixParameter(4,19828.9);
    Mean->FixParameter(5,1);
    Mean->FixParameter(6,1.8756129425);

  }

  double mean = Mean->Eval(p);
  Mean->Delete();

  double Resolution = 0.0;

  if(((isProton == true) || (isAntiProton == true))	&& (MetaLHC16 == true)) Resolution = 0.126;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && (MetaLHC16 == true)) Resolution = 9.14588e-02;

  if(((isProton == true) || (isAntiProton == true))	&& (MetaLHC17 == true)) Resolution = 1.34216e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && (MetaLHC17 == true)) Resolution = 9.00246e-02;

  if(((isProton == true) || (isAntiProton == true))	&& (MetaLHC18 == true)) Resolution = 1.32506e-01;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && (MetaLHC18 == true)) Resolution = 8.82121e-02;

  if(((isProton == true) || (isAntiProton == true))	&& ((LHC18q == true) || (LHC18r == true))) Resolution = 0.10;
  if(((isDeuteron == true) || (isAntiDeuteron == true)) && ((LHC18q == true) || (LHC18r == true))) Resolution = 0.10;

  double ScaleFactor = 1.0-(Resolution);
  double sigma = (mean*ScaleFactor) - mean;

  SigmaParticle = (mean - SignalITS) / (sigma);

  return SigmaParticle;

} // end of CalculateSigmadEdxITS



