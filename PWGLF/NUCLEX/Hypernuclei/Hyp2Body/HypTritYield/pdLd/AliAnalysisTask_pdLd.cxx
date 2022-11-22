#include "TChain.h" 
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
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
#include "AliEventPoolManager.h"

#include "AliKFVertex.h"
#include "AliKFParticleBase.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <string>

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenHijingEventHeader.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TRandom3.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTask_pdLd.h"

using namespace std;
ClassImp(AliAnalysisTask_pdLd) 
ClassImp(AliAODTrackTiny) 






AliAnalysisTask_pdLd::AliAnalysisTask_pdLd() : AliAnalysisTaskSE(),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fEventPoolManager(0),
  fStack(0),
  fHistList_Event(0),
  fHistList_Proton(0),
  fHistList_Deuteron(0),
  fHistList_ProtonDeuteron(0),
  fHistList_AntiProton(0),
  fHistList_AntiDeuteron(0),
  fHistList_AntiProtonAntiDeuteron(0),
  fCollisionSystem(0),
  fHist_Event_CutCounter(0),
  fHist_Event_PrimVertexZ(0),
  fHist_Event_Centrality(0),
  fHist_Proton_CutCounter(0),
  fHist_Proton_pT(0),
  fHist_Proton_p(0),
  fHist_Proton_pTPC(0),
  fHist_Proton_Eta(0),
  fHist_Proton_Phi(0),
  fHist_Proton_TPC_nCluster(0),
  fHist_Proton_TPC_CrossedRows(0),
  fHist_Proton_TPC_RatioRowsCluster(0),
  fHist_Proton_TPC_SharedCluster(0),
  fHist_Proton_TPC_Chi2overNDF(0),
  fHist_Proton_ITS_nCluster(0),
  fHist_Proton_TOF_Signal(0),
  fHist_Proton_DCAxy(0),
  fHist_Proton_DCAz(0),
  fHist_Proton_TPC_nSigma_pT(0),
  fHist_Proton_TPC_nSigma_p(0),
  fHist_Proton_TOF_nSigma_pT(0),
  fHist_Proton_TOF_nSigma_p(0),
  fHist_Proton_ITS_nSigma_pT(0),
  fHist_Proton_ITS_nSigma_p(0),
  fHist_Proton_TPC_dEdx_pT(0),
  fHist_Proton_TPC_dEdx_p(0),
  fHist_Proton_TOF_Beta_pT(0),
  fHist_Proton_TOF_Beta_p(0),
  fHist_Proton_TOF_MassSquare_pT(0),
  fHist_Proton_TOF_MassSquare_p(0),
  fHist_Proton_ITS_dEdx_pT(0),
  fHist_Proton_ITS_dEdx_p(0),
  fHist_Deuteron_CutCounter(0),
  fHist_Deuteron_pT(0),
  fHist_Deuteron_p(0),
  fHist_Deuteron_pTPC(0),
  fHist_Deuteron_Eta(0),
  fHist_Deuteron_Phi(0),
  fHist_Deuteron_TPC_nCluster(0),
  fHist_Deuteron_TPC_CrossedRows(0),
  fHist_Deuteron_TPC_RatioRowsCluster(0),
  fHist_Deuteron_TPC_SharedCluster(0),
  fHist_Deuteron_TPC_Chi2overNDF(0),
  fHist_Deuteron_ITS_nCluster(0),
  fHist_Deuteron_TOF_Signal(0),
  fHist_Deuteron_DCAxy(0),
  fHist_Deuteron_DCAz(0),
  fHist_Deuteron_TPC_nSigma_pT(0),
  fHist_Deuteron_TPC_nSigma_p(0),
  fHist_Deuteron_TOF_nSigma_pT(0),
  fHist_Deuteron_TOF_nSigma_p(0),
  fHist_Deuteron_TOF_nSigmaMassSq_pT(0),
  fHist_Deuteron_TOF_nSigmaMassSq_p(0),
  fHist_Deuteron_ITS_nSigma_pT(0),
  fHist_Deuteron_ITS_nSigma_p(0),
  fHist_Deuteron_TPC_dEdx_pT(0),
  fHist_Deuteron_TPC_dEdx_p(0),
  fHist_Deuteron_TOF_Beta_pT(0),
  fHist_Deuteron_TOF_Beta_p(0),
  fHist_Deuteron_TOF_MassSquare_pT(0),
  fHist_Deuteron_TOF_MassSquare_p(0),
  fHist_Deuteron_ITS_dEdx_pT(0),
  fHist_Deuteron_ITS_dEdx_p(0),
  fHist_ProtonDeuteron_CutCounter(0),
  fHist_ProtonDeuteron_SED(0),
  fHist_ProtonDeuteron_MED(0),
  fHist_ProtonDeuteron_RPD(0),
  fHist_ProtonDeuteron_PairsPerEvent(0),
  fHist_ProtonDeuteron_EventsForMixing(0),
  fHist_ProtonDeuteron_AngleOfPairs(0),
  fHist_ProtonDeuteron_PairMultiplicity(0),
  fHist_ProtonDeuteron_pT(0),
  fHist_ProtonDeuteron_Eta(0),
  fHist_ProtonDeuteron_Centrality(0),
  fHist_ProtonDeuteron_VertexZ(0),
  fHist_ProtonDeuteron_UsedEventsInPool(0),
  fHist_AntiProton_CutCounter(0),
  fHist_AntiProton_pT(0),
  fHist_AntiProton_p(0),
  fHist_AntiProton_pTPC(0),
  fHist_AntiProton_Eta(0),
  fHist_AntiProton_Phi(0),
  fHist_AntiProton_TPC_nCluster(0),
  fHist_AntiProton_TPC_CrossedRows(0),
  fHist_AntiProton_TPC_RatioRowsCluster(0),
  fHist_AntiProton_TPC_SharedCluster(0),
  fHist_AntiProton_TPC_Chi2overNDF(0),
  fHist_AntiProton_ITS_nCluster(0),
  fHist_AntiProton_TOF_Signal(0),
  fHist_AntiProton_DCAxy(0),
  fHist_AntiProton_DCAz(0),
  fHist_AntiProton_TPC_nSigma_pT(0),
  fHist_AntiProton_TPC_nSigma_p(0),
  fHist_AntiProton_TOF_nSigma_pT(0),
  fHist_AntiProton_TOF_nSigma_p(0),
  fHist_AntiProton_ITS_nSigma_pT(0),
  fHist_AntiProton_ITS_nSigma_p(0),
  fHist_AntiProton_TPC_dEdx_pT(0),
  fHist_AntiProton_TPC_dEdx_p(0),
  fHist_AntiProton_TOF_Beta_pT(0),
  fHist_AntiProton_TOF_Beta_p(0),
  fHist_AntiProton_TOF_MassSquare_pT(0),
  fHist_AntiProton_TOF_MassSquare_p(0),
  fHist_AntiProton_ITS_dEdx_pT(0),
  fHist_AntiProton_ITS_dEdx_p(0),
  fHist_AntiDeuteron_CutCounter(0),
  fHist_AntiDeuteron_pT(0),
  fHist_AntiDeuteron_p(0),
  fHist_AntiDeuteron_pTPC(0),
  fHist_AntiDeuteron_Eta(0),
  fHist_AntiDeuteron_Phi(0),
  fHist_AntiDeuteron_TPC_nCluster(0),
  fHist_AntiDeuteron_TPC_CrossedRows(0),
  fHist_AntiDeuteron_TPC_RatioRowsCluster(0),
  fHist_AntiDeuteron_TPC_SharedCluster(0),
  fHist_AntiDeuteron_TPC_Chi2overNDF(0),
  fHist_AntiDeuteron_ITS_nCluster(0),
  fHist_AntiDeuteron_TOF_Signal(0),
  fHist_AntiDeuteron_DCAxy(0),
  fHist_AntiDeuteron_DCAz(0),
  fHist_AntiDeuteron_TPC_nSigma_pT(0),
  fHist_AntiDeuteron_TPC_nSigma_p(0),
  fHist_AntiDeuteron_TOF_nSigma_pT(0),
  fHist_AntiDeuteron_TOF_nSigma_p(0),
  fHist_AntiDeuteron_TOF_nSigmaMassSq_pT(0),
  fHist_AntiDeuteron_TOF_nSigmaMassSq_p(0),
  fHist_AntiDeuteron_ITS_nSigma_pT(0),
  fHist_AntiDeuteron_ITS_nSigma_p(0),
  fHist_AntiDeuteron_TPC_dEdx_pT(0),
  fHist_AntiDeuteron_TPC_dEdx_p(0),
  fHist_AntiDeuteron_TOF_Beta_pT(0),
  fHist_AntiDeuteron_TOF_Beta_p(0),
  fHist_AntiDeuteron_TOF_MassSquare_pT(0),
  fHist_AntiDeuteron_TOF_MassSquare_p(0),
  fHist_AntiDeuteron_ITS_dEdx_pT(0),
  fHist_AntiDeuteron_ITS_dEdx_p(0),
  fHist_AntiProtonAntiDeuteron_CutCounter(0),
  fHist_AntiProtonAntiDeuteron_SED(0),
  fHist_AntiProtonAntiDeuteron_MED(0),
  fHist_AntiProtonAntiDeuteron_RPD(0),
  fHist_AntiProtonAntiDeuteron_PairsPerEvent(0),
  fHist_AntiProtonAntiDeuteron_EventsForMixing(0),
  fHist_AntiProtonAntiDeuteron_AngleOfPairs(0),
  fHist_AntiProtonAntiDeuteron_PairMultiplicity(0),
  fHist_AntiProtonAntiDeuteron_pT(0),
  fHist_AntiProtonAntiDeuteron_Eta(0),
  fHist_AntiProtonAntiDeuteron_Centrality(0),
  fHist_AntiProtonAntiDeuteron_VertexZ(0),
  fHist_AntiProtonAntiDeuteron_UsedEventsInPool(0),
  ProtonTrackArray(0),
  DeuteronTrackArray(0),
  Lambdav0Array(0),
  AntiProtonTrackArray(0),
  AntiDeuteronTrackArray(0),
  AntiLambdav0Array(0)
{


}



AliAnalysisTask_pdLd::AliAnalysisTask_pdLd(const char *name,int CollisionSystem) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fEventPoolManager(0),
  fStack(0),
  fHistList_Event(0),
  fHistList_Proton(0),
  fHistList_Deuteron(0),
  fHistList_ProtonDeuteron(0),
  fHistList_AntiProton(0),
  fHistList_AntiDeuteron(0),
  fHistList_AntiProtonAntiDeuteron(0),
  fCollisionSystem(CollisionSystem),
  fHist_Event_CutCounter(0),
  fHist_Event_PrimVertexZ(0),
  fHist_Event_Centrality(0),
  fHist_Proton_CutCounter(0),
  fHist_Proton_pT(0),
  fHist_Proton_p(0),
  fHist_Proton_pTPC(0),
  fHist_Proton_Eta(0),
  fHist_Proton_Phi(0),
  fHist_Proton_TPC_nCluster(0),
  fHist_Proton_TPC_CrossedRows(0),
  fHist_Proton_TPC_RatioRowsCluster(0),
  fHist_Proton_TPC_SharedCluster(0),
  fHist_Proton_TPC_Chi2overNDF(0),
  fHist_Proton_ITS_nCluster(0),
  fHist_Proton_TOF_Signal(0),
  fHist_Proton_DCAxy(0),
  fHist_Proton_DCAz(0),
  fHist_Proton_TPC_nSigma_pT(0),
  fHist_Proton_TPC_nSigma_p(0),
  fHist_Proton_TOF_nSigma_pT(0),
  fHist_Proton_TOF_nSigma_p(0),
  fHist_Proton_ITS_nSigma_pT(0),
  fHist_Proton_ITS_nSigma_p(0),
  fHist_Proton_TPC_dEdx_pT(0),
  fHist_Proton_TPC_dEdx_p(0),
  fHist_Proton_TOF_Beta_pT(0),
  fHist_Proton_TOF_Beta_p(0),
  fHist_Proton_TOF_MassSquare_pT(0),
  fHist_Proton_TOF_MassSquare_p(0),
  fHist_Proton_ITS_dEdx_pT(0),
  fHist_Proton_ITS_dEdx_p(0),
  fHist_Deuteron_CutCounter(0),
  fHist_Deuteron_pT(0),
  fHist_Deuteron_p(0),
  fHist_Deuteron_pTPC(0),
  fHist_Deuteron_Eta(0),
  fHist_Deuteron_Phi(0),
  fHist_Deuteron_TPC_nCluster(0),
  fHist_Deuteron_TPC_CrossedRows(0),
  fHist_Deuteron_TPC_RatioRowsCluster(0),
  fHist_Deuteron_TPC_SharedCluster(0),
  fHist_Deuteron_TPC_Chi2overNDF(0),
  fHist_Deuteron_ITS_nCluster(0),
  fHist_Deuteron_TOF_Signal(0),
  fHist_Deuteron_DCAxy(0),
  fHist_Deuteron_DCAz(0),
  fHist_Deuteron_TPC_nSigma_pT(0),
  fHist_Deuteron_TPC_nSigma_p(0),
  fHist_Deuteron_TOF_nSigma_pT(0),
  fHist_Deuteron_TOF_nSigma_p(0),
  fHist_Deuteron_TOF_nSigmaMassSq_pT(0),
  fHist_Deuteron_TOF_nSigmaMassSq_p(0),
  fHist_Deuteron_ITS_nSigma_pT(0),
  fHist_Deuteron_ITS_nSigma_p(0),
  fHist_Deuteron_TPC_dEdx_pT(0),
  fHist_Deuteron_TPC_dEdx_p(0),
  fHist_Deuteron_TOF_Beta_pT(0),
  fHist_Deuteron_TOF_Beta_p(0),
  fHist_Deuteron_TOF_MassSquare_pT(0),
  fHist_Deuteron_TOF_MassSquare_p(0),
  fHist_Deuteron_ITS_dEdx_pT(0),
  fHist_Deuteron_ITS_dEdx_p(0),
  fHist_ProtonDeuteron_CutCounter(0),
  fHist_ProtonDeuteron_SED(0),
  fHist_ProtonDeuteron_MED(0),
  fHist_ProtonDeuteron_RPD(0),
  fHist_ProtonDeuteron_PairsPerEvent(0),
  fHist_ProtonDeuteron_EventsForMixing(0),
  fHist_ProtonDeuteron_AngleOfPairs(0),
  fHist_ProtonDeuteron_PairMultiplicity(0),
  fHist_ProtonDeuteron_pT(0),
  fHist_ProtonDeuteron_Eta(0),
  fHist_ProtonDeuteron_Centrality(0),
  fHist_ProtonDeuteron_VertexZ(0),
  fHist_ProtonDeuteron_UsedEventsInPool(0),
  fHist_AntiProton_CutCounter(0),
  fHist_AntiProton_pT(0),
  fHist_AntiProton_p(0),
  fHist_AntiProton_pTPC(0),
  fHist_AntiProton_Eta(0),
  fHist_AntiProton_Phi(0),
  fHist_AntiProton_TPC_nCluster(0),
  fHist_AntiProton_TPC_CrossedRows(0),
  fHist_AntiProton_TPC_RatioRowsCluster(0),
  fHist_AntiProton_TPC_SharedCluster(0),
  fHist_AntiProton_TPC_Chi2overNDF(0),
  fHist_AntiProton_ITS_nCluster(0),
  fHist_AntiProton_TOF_Signal(0),
  fHist_AntiProton_DCAxy(0),
  fHist_AntiProton_DCAz(0),
  fHist_AntiProton_TPC_nSigma_pT(0),
  fHist_AntiProton_TPC_nSigma_p(0),
  fHist_AntiProton_TOF_nSigma_pT(0),
  fHist_AntiProton_TOF_nSigma_p(0),
  fHist_AntiProton_ITS_nSigma_pT(0),
  fHist_AntiProton_ITS_nSigma_p(0),
  fHist_AntiProton_TPC_dEdx_pT(0),
  fHist_AntiProton_TPC_dEdx_p(0),
  fHist_AntiProton_TOF_Beta_pT(0),
  fHist_AntiProton_TOF_Beta_p(0),
  fHist_AntiProton_TOF_MassSquare_pT(0),
  fHist_AntiProton_TOF_MassSquare_p(0),
  fHist_AntiProton_ITS_dEdx_pT(0),
  fHist_AntiProton_ITS_dEdx_p(0),
  fHist_AntiDeuteron_CutCounter(0),
  fHist_AntiDeuteron_pT(0),
  fHist_AntiDeuteron_p(0),
  fHist_AntiDeuteron_pTPC(0),
  fHist_AntiDeuteron_Eta(0),
  fHist_AntiDeuteron_Phi(0),
  fHist_AntiDeuteron_TPC_nCluster(0),
  fHist_AntiDeuteron_TPC_CrossedRows(0),
  fHist_AntiDeuteron_TPC_RatioRowsCluster(0),
  fHist_AntiDeuteron_TPC_SharedCluster(0),
  fHist_AntiDeuteron_TPC_Chi2overNDF(0),
  fHist_AntiDeuteron_ITS_nCluster(0),
  fHist_AntiDeuteron_TOF_Signal(0),
  fHist_AntiDeuteron_DCAxy(0),
  fHist_AntiDeuteron_DCAz(0),
  fHist_AntiDeuteron_TPC_nSigma_pT(0),
  fHist_AntiDeuteron_TPC_nSigma_p(0),
  fHist_AntiDeuteron_TOF_nSigma_pT(0),
  fHist_AntiDeuteron_TOF_nSigma_p(0),
  fHist_AntiDeuteron_TOF_nSigmaMassSq_pT(0),
  fHist_AntiDeuteron_TOF_nSigmaMassSq_p(0),
  fHist_AntiDeuteron_ITS_nSigma_pT(0),
  fHist_AntiDeuteron_ITS_nSigma_p(0),
  fHist_AntiDeuteron_TPC_dEdx_pT(0),
  fHist_AntiDeuteron_TPC_dEdx_p(0),
  fHist_AntiDeuteron_TOF_Beta_pT(0),
  fHist_AntiDeuteron_TOF_Beta_p(0),
  fHist_AntiDeuteron_TOF_MassSquare_pT(0),
  fHist_AntiDeuteron_TOF_MassSquare_p(0),
  fHist_AntiDeuteron_ITS_dEdx_pT(0),
  fHist_AntiDeuteron_ITS_dEdx_p(0),
  fHist_AntiProtonAntiDeuteron_CutCounter(0),
  fHist_AntiProtonAntiDeuteron_SED(0),
  fHist_AntiProtonAntiDeuteron_MED(0),
  fHist_AntiProtonAntiDeuteron_RPD(0),
  fHist_AntiProtonAntiDeuteron_PairsPerEvent(0),
  fHist_AntiProtonAntiDeuteron_EventsForMixing(0),
  fHist_AntiProtonAntiDeuteron_AngleOfPairs(0),
  fHist_AntiProtonAntiDeuteron_PairMultiplicity(0),
  fHist_AntiProtonAntiDeuteron_pT(0),
  fHist_AntiProtonAntiDeuteron_Eta(0),
  fHist_AntiProtonAntiDeuteron_Centrality(0),
  fHist_AntiProtonAntiDeuteron_VertexZ(0),
  fHist_AntiProtonAntiDeuteron_UsedEventsInPool(0),
  ProtonTrackArray(0),
  DeuteronTrackArray(0),
  Lambdav0Array(0),
  AntiProtonTrackArray(0),
  AntiDeuteronTrackArray(0),
  AntiLambdav0Array(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TList::Class());
  DefineOutput(5,TList::Class());
  DefineOutput(6,TList::Class());
  DefineOutput(7,TList::Class());

}

  
AliAnalysisTask_pdLd::~AliAnalysisTask_pdLd()
{

  if(fHistList_Event)
    {
      fHistList_Event->Clear();
      delete fHistList_Event;
    }

  if(fHistList_Proton)
    {
      fHistList_Proton->Clear();
      delete fHistList_Proton;
    }

  if(fHistList_Deuteron)
    {
      fHistList_Deuteron->Clear();
      delete fHistList_Deuteron;
    }

  if(fHistList_ProtonDeuteron)
    {
      fHistList_ProtonDeuteron->Clear();
      delete fHistList_ProtonDeuteron;
    }

  if(fHistList_AntiProton)
    {
      fHistList_AntiProton->Clear();
      delete fHistList_AntiProton;
    }

  if(fHistList_AntiDeuteron)
    {
      fHistList_AntiDeuteron->Clear();
      delete fHistList_AntiDeuteron;
    }

  if(fHistList_AntiProtonAntiDeuteron)
    {
      fHistList_AntiProtonAntiDeuteron->Clear();
      delete fHistList_AntiProtonAntiDeuteron;
    }


}







void AliAnalysisTask_pdLd::UserCreateOutputObjects()
{

  fHistList_Event = new TList();
  fHistList_Event->SetOwner();

  fHistList_Proton = new TList();
  fHistList_Proton->SetOwner();

  fHistList_Deuteron = new TList();
  fHistList_Deuteron->SetOwner();

  fHistList_ProtonDeuteron = new TList();
  fHistList_ProtonDeuteron->SetOwner();

  fHistList_AntiProton = new TList();
  fHistList_AntiProton->SetOwner();

  fHistList_AntiDeuteron = new TList();
  fHistList_AntiDeuteron->SetOwner();

  fHistList_AntiProtonAntiDeuteron = new TList();
  fHistList_AntiProtonAntiDeuteron->SetOwner();


  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ PoolManager for event-mixing ++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  const int nCentralityBins = 20; // number of centrality bins of the centrality array
  const int nZvtxBins	    = 20; // number of z-vertex bins of the z-vertex array

  double CentralityBins[nCentralityBins+1];

  double CentralityBins_Central[nCentralityBins+1] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};
  double CentralityBins_SemiCentral[nCentralityBins+1] = {30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0};

  if(fCollisionSystem == 1) for(int i = 0; i <= nCentralityBins; i++){CentralityBins[i] = CentralityBins_Central[i];}
  if(fCollisionSystem == 2) for(int i = 0; i <= nCentralityBins; i++){CentralityBins[i] = CentralityBins_SemiCentral[i];}


  double ZvtxBins[nZvtxBins+1] = {-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};


  int PoolSize = 200;	    // maximum number of events in the pool (-1 means no limit)
			    // int nPools = nCentralityBins * nZvtxBins;
			    // int nEventsMax = nPools * PoolSize;

  int TrackDepth = 5000;   // maximum number of tracks in one bin of the pool?

  fEventPoolManager = new AliEventPoolManager(PoolSize,TrackDepth,nCentralityBins,CentralityBins,nZvtxBins,ZvtxBins);
  fEventPoolManager->Validate();

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for events +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_Event_CutCounter = new TH1F("fHist_Event_CutCounter","number of events containing selected...",7,0.0,7.0);
  fHist_Event_CutCounter->GetYaxis()->SetTitle("number of events");
  fHist_Event_CutCounter->GetXaxis()->SetNdivisions(107,false);
  fHist_Event_CutCounter->GetXaxis()->CenterLabels();
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(1,"total");
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(2,"p");
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(3,"d");
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(4,"p-d");
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(5,"#bar{p}");
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(6,"#bar{d}");
  fHist_Event_CutCounter->GetXaxis()->SetBinLabel(7,"#bar{p}-#bar{d}");
  fHist_Event_CutCounter->SetStats(0);
  fHistList_Event->Add(fHist_Event_CutCounter);

  fHist_Event_PrimVertexZ = new TH1F("fHist_Event_PrimVertexZ","z-position of primary vertex",240,-12.0,+12.0);
  fHist_Event_PrimVertexZ->GetXaxis()->SetTitle("z-position of primary vertex (cm)");
  fHist_Event_PrimVertexZ->GetYaxis()->SetTitle("counts");
  fHistList_Event->Add(fHist_Event_PrimVertexZ);

  fHist_Event_Centrality = new TH1F("fHist_Event_Centrality","event centrality",100,0.0,100.0);
  fHist_Event_Centrality->GetXaxis()->SetTitle("centrality (%)");
  fHist_Event_Centrality->GetYaxis()->SetTitle("counts");
  fHistList_Event->Add(fHist_Event_Centrality);
  
  
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for protons ++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_Proton_CutCounter = new TH1F("fHist_Proton_CutCounter","cut counter of proton selection",20,0.0,20.0);
  fHist_Proton_CutCounter->GetYaxis()->SetTitle("counts");
  fHist_Proton_CutCounter->GetXaxis()->SetNdivisions(120,false);
  fHist_Proton_CutCounter->GetXaxis()->CenterLabels();
  fHist_Proton_CutCounter->GetXaxis()->SetLabelSize(0.025);
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(1,"no cuts");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(2,"TPC PID status");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(3,"n#sigma_{TPC}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(4,"DCA_{xy}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(5,"DCA_{z}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(6,"TOF PID status");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(7,"n#sigma_{TOF}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(8,"reject smaller #sigma");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(9,"FilterBit");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(10,"#it{p}_{T}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(11,"charge");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(12,"#eta");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(13,"n_{Cluster (TPC)}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(14,"n_{CrossedRows (TPC)}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(15,"n_{SharedCluster (TPC)}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(16,"#frac{n_{CrossedRows}}{n_{FindableCluster}}(TPC)");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(17,"low-#it{p}_{T} #pi");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(18,"ITS PID status");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(19,"n#sigma_{ITS}");
  fHist_Proton_CutCounter->GetXaxis()->SetBinLabel(20,"n_{Cluster} (ITS)");
  fHist_Proton_CutCounter->SetStats(0);
  fHistList_Proton->Add(fHist_Proton_CutCounter);

  fHist_Proton_pT = new TH1F("fHist_Proton_pT","#it{p}_{T} distribution of protons",240,0.0,6.0);
  fHist_Proton_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_pT->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_pT);

  fHist_Proton_p = new TH1F("fHist_Proton_p","#it{p} distribution of protons",240,0.0,6.0);
  fHist_Proton_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_Proton_p->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_p);

  fHist_Proton_pTPC = new TH1F("fHist_Proton_pTPC","#it{p}_{TPC} distribution of protons",240,0.0,6.0);
  fHist_Proton_pTPC->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Proton_pTPC->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_pTPC);

  fHist_Proton_Eta = new TH1F("fHist_Proton_Eta","#it{#eta} distribution of protons",200,-1.5,+1.5);
  fHist_Proton_Eta->GetXaxis()->SetTitle("#it{#eta}");
  fHist_Proton_Eta->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_Eta);

  fHist_Proton_Phi = new TH1F("fHist_Proton_Phi","#it{#phi} distribution of protons",360,0.0,360.0);
  fHist_Proton_Phi->GetXaxis()->SetTitle("#it{#phi} (#circ)");
  fHist_Proton_Phi->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_Phi);

  fHist_Proton_TPC_nCluster = new TH1F("fHist_Proton_TPC_nCluster","Number of TPC clusters of proton tracks",160,0.0,160.0);
  fHist_Proton_TPC_nCluster->GetXaxis()->SetTitle("n_{Cluster,TPC}");
  fHist_Proton_TPC_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_TPC_nCluster);

  fHist_Proton_TPC_CrossedRows = new TH1F("fHist_Proton_TPC_CrossedRows","Number of crossed TPC rows of proton tracks",160,0.0,160.0);
  fHist_Proton_TPC_CrossedRows->GetXaxis()->SetTitle("n_{Rows,TPC}");
  fHist_Proton_TPC_CrossedRows->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_TPC_CrossedRows);

  fHist_Proton_TPC_RatioRowsCluster = new TH1F("fHist_Proton_TPC_RatioRowsCluster","Crossed TPC rows / Findable clusters for protons",150,0.0,1.5);
  fHist_Proton_TPC_RatioRowsCluster->GetXaxis()->SetTitle("n_{Rows,TPC}/n_{Cluster,TPC}");
  fHist_Proton_TPC_RatioRowsCluster->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_TPC_RatioRowsCluster);

  fHist_Proton_TPC_SharedCluster = new TH1F("fHist_Proton_TPC_SharedCluster","Number of shared TPC clusters of proton tracks",160,0.0,160.0);
  fHist_Proton_TPC_SharedCluster->GetXaxis()->SetTitle("n_{Shared Clusters,TPC}");
  fHist_Proton_TPC_SharedCluster->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_TPC_SharedCluster);

  fHist_Proton_TPC_Chi2overNDF = new TH1F("fHist_Proton_TPC_Chi2overNDF","#it{#chi}^{2}/ndf of TPC proton tracks",100,0.0,5.0);
  fHist_Proton_TPC_Chi2overNDF->GetXaxis()->SetTitle("#it{#chi}^{2}/ndf");
  fHist_Proton_TPC_Chi2overNDF->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_TPC_Chi2overNDF);

  fHist_Proton_ITS_nCluster = new TH1F("fHist_Proton_ITS_nCluster","Number of ITS clusters of proton tracks",20,0.0,20.0);
  fHist_Proton_ITS_nCluster->GetXaxis()->SetTitle("n_{Cluster} ITS");
  fHist_Proton_ITS_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_ITS_nCluster);

  fHist_Proton_TOF_Signal = new TH1F("fHist_Proton_TOF_Signal","TOF signal of proton tracks",200,0.0,200.0);
  fHist_Proton_TOF_Signal->GetXaxis()->SetTitle("TOF signal (ns)");
  fHist_Proton_TOF_Signal->GetYaxis()->SetTitle("counts");
  fHistList_Proton->Add(fHist_Proton_TOF_Signal);

  fHist_Proton_DCAxy = new TH2F("fHist_Proton_DCAxy","DCA_{xy} distribution of proton tracks",240,0.0,6.0,1000,-3.0,+3.0);
  fHist_Proton_DCAxy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_DCAxy->GetYaxis()->SetTitle("DCA_{xy} (cm)");
  fHistList_Proton->Add(fHist_Proton_DCAxy);

  fHist_Proton_DCAz = new TH2F("fHist_Proton_DCAz","DCA_{z} distribution of proton tracks",240,0.0,6.0,1000,-3.0,+3.0);
  fHist_Proton_DCAz->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_DCAz->GetYaxis()->SetTitle("DCA_{z} (cm)");
  fHistList_Proton->Add(fHist_Proton_DCAz);

  fHist_Proton_TPC_nSigma_pT = new TH2F("fHist_Proton_TPC_nSigma_pT","n#it{#sigma}_{TPC} of proton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Proton_TPC_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_TPC_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_Proton->Add(fHist_Proton_TPC_nSigma_pT);

  fHist_Proton_TPC_nSigma_p = new TH2F("fHist_Proton_TPC_nSigma_p","n#it{#sigma}_{TPC} of proton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Proton_TPC_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Proton_TPC_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_Proton->Add(fHist_Proton_TPC_nSigma_p);

  fHist_Proton_TOF_nSigma_pT = new TH2F("fHist_Proton_TOF_nSigma_pT","n#it{#sigma}_{TOF} of proton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Proton_TOF_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_TOF_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_Proton->Add(fHist_Proton_TOF_nSigma_pT);

  fHist_Proton_TOF_nSigma_p = new TH2F("fHist_Proton_TOF_nSigma_p","n#it{#sigma}_{TOF} of proton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Proton_TOF_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Proton_TOF_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_Proton->Add(fHist_Proton_TOF_nSigma_p);

  fHist_Proton_ITS_nSigma_pT = new TH2F("fHist_Proton_ITS_nSigma_pT","n#it{#sigma}_{ITS} of proton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Proton_ITS_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_ITS_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_Proton->Add(fHist_Proton_ITS_nSigma_pT);

  fHist_Proton_ITS_nSigma_p = new TH2F("fHist_Proton_ITS_nSigma_p","n#it{#sigma}_{ITS} of proton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Proton_ITS_nSigma_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_Proton_ITS_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_Proton->Add(fHist_Proton_ITS_nSigma_p);

  fHist_Proton_TPC_dEdx_pT = new TH2F("fHist_Proton_TPC_dEdx_pT","d#it{E}/d#it{x} of protons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_Proton_TPC_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_TPC_dEdx_pT->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_Proton->Add(fHist_Proton_TPC_dEdx_pT);

  fHist_Proton_TPC_dEdx_p = new TH2F("fHist_Proton_TPC_dEdx_p","d#it{E}/d#it{x} of protons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_Proton_TPC_dEdx_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Proton_TPC_dEdx_p->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_Proton->Add(fHist_Proton_TPC_dEdx_p);

  fHist_Proton_TOF_Beta_pT = new TH2F("fHist_Proton_TOF_Beta_pT","#it{#beta} of protons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_Proton_TOF_Beta_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_TOF_Beta_pT->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_Proton->Add(fHist_Proton_TOF_Beta_pT);

  fHist_Proton_TOF_Beta_p = new TH2F("fHist_Proton_TOF_Beta_p","#it{#beta} of protons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_Proton_TOF_Beta_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Proton_TOF_Beta_p->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_Proton->Add(fHist_Proton_TOF_Beta_p);

  fHist_Proton_TOF_MassSquare_pT = new TH2F("fHist_Proton_TOF_MassSquare_pT","#it{m}^{2} of protons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_Proton_TOF_MassSquare_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_TOF_MassSquare_pT->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_Proton->Add(fHist_Proton_TOF_MassSquare_pT);

  fHist_Proton_TOF_MassSquare_p = new TH2F("fHist_Proton_TOF_MassSquare_p","#it{m}^{2} of protons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_Proton_TOF_MassSquare_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Proton_TOF_MassSquare_p->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_Proton->Add(fHist_Proton_TOF_MassSquare_p);

  fHist_Proton_ITS_dEdx_pT = new TH2F("fHist_Proton_ITS_dEdx_pT","d#it{E}/d#it{x} of protons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_Proton_ITS_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Proton_ITS_dEdx_pT->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_Proton->Add(fHist_Proton_ITS_dEdx_pT);

  fHist_Proton_ITS_dEdx_p = new TH2F("fHist_Proton_ITS_dEdx_p","d#it{E}/d#it{x} of protons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_Proton_ITS_dEdx_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_Proton_ITS_dEdx_p->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_Proton->Add(fHist_Proton_ITS_dEdx_p);






  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for deuterons ++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_Deuteron_CutCounter = new TH1F("fHist_Deuteron_CutCounter","cut counter of deuteron selection",21,0.0,21.0);
  fHist_Deuteron_CutCounter->GetYaxis()->SetTitle("counts");
  fHist_Deuteron_CutCounter->GetXaxis()->SetNdivisions(121,false);
  fHist_Deuteron_CutCounter->GetXaxis()->CenterLabels();
  fHist_Deuteron_CutCounter->GetXaxis()->SetLabelSize(0.025);
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(1,"no cuts");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(2,"TPC PID status");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(3,"n#sigma_{TPC}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(4,"DCA_{xy}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(5,"DCA_{z}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(6,"TOF PID status");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(7,"n#sigma_{TOF}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(8,"reject smaller #sigma");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(9,"n#sigma_{TOF} m^{2}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(10,"FilterBit");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(11,"#it{p}_{T}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(12,"charge");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(13,"#eta");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(14,"n_{Cluster (TPC)}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(15,"n_{CrossedRows (TPC)}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(16,"n_{SharedCluster (TPC)}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(17,"#frac{n_{CrossedRows}}{n_{FindableCluster}}(TPC)");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(18,"low-#it{p}_{T} #pi");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(19,"ITS PID status");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(20,"n#sigma_{ITS}");
  fHist_Deuteron_CutCounter->GetXaxis()->SetBinLabel(21,"n_{Cluster} (ITS)");
  fHist_Deuteron_CutCounter->SetStats(0);
  fHistList_Deuteron->Add(fHist_Deuteron_CutCounter);

  fHist_Deuteron_pT = new TH1F("fHist_Deuteron_pT","#it{p}_{T} distribution of deuterons",240,0.0,6.0);
  fHist_Deuteron_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_pT->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_pT);

  fHist_Deuteron_p = new TH1F("fHist_Deuteron_p","#it{p} distribution of deuterons",240,0.0,6.0);
  fHist_Deuteron_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_Deuteron_p->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_p);

  fHist_Deuteron_pTPC = new TH1F("fHist_Deuteron_pTPC","#it{p}_{TPC} distribution of deuterons",240,0.0,6.0);
  fHist_Deuteron_pTPC->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_pTPC->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_pTPC);

  fHist_Deuteron_Eta = new TH1F("fHist_Deuteron_Eta","#it{#eta} distribution of deuterons",200,-1.5,+1.5);
  fHist_Deuteron_Eta->GetXaxis()->SetTitle("#it{#eta}");
  fHist_Deuteron_Eta->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_Eta);

  fHist_Deuteron_Phi = new TH1F("fHist_Deuteron_Phi","#it{#phi} distribution of deuterons",360,0.0,360.0);
  fHist_Deuteron_Phi->GetXaxis()->SetTitle("#it{#phi} (#circ)");
  fHist_Deuteron_Phi->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_Phi);

  fHist_Deuteron_TPC_nCluster = new TH1F("fHist_Deuteron_TPC_nCluster","Number of TPC clusters of deuteron tracks",160,0.0,160.0);
  fHist_Deuteron_TPC_nCluster->GetXaxis()->SetTitle("n_{Cluster,TPC}");
  fHist_Deuteron_TPC_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_nCluster);

  fHist_Deuteron_TPC_CrossedRows = new TH1F("fHist_Deuteron_TPC_CrossedRows","Number of crossed TPC rows of deuteron tracks",160,0.0,160.0);
  fHist_Deuteron_TPC_CrossedRows->GetXaxis()->SetTitle("n_{Rows,TPC}");
  fHist_Deuteron_TPC_CrossedRows->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_CrossedRows);

  fHist_Deuteron_TPC_RatioRowsCluster = new TH1F("fHist_Deuteron_TPC_RatioRowsCluster","Crossed TPC rows / Findable clusters for deuterons",150,0.0,1.5);
  fHist_Deuteron_TPC_RatioRowsCluster->GetXaxis()->SetTitle("n_{Rows,TPC}/n_{Cluster,TPC}");
  fHist_Deuteron_TPC_RatioRowsCluster->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_RatioRowsCluster);

  fHist_Deuteron_TPC_SharedCluster = new TH1F("fHist_Deuteron_TPC_SharedCluster","Number of shared TPC clusters of deuteron tracks",160,0.0,160.0);
  fHist_Deuteron_TPC_SharedCluster->GetXaxis()->SetTitle("n_{Shared Clusters,TPC}");
  fHist_Deuteron_TPC_SharedCluster->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_SharedCluster);

  fHist_Deuteron_TPC_Chi2overNDF = new TH1F("fHist_Deuteron_TPC_Chi2overNDF","#it{#chi}^{2}/ndf of TPC deuteron tracks",100,0.0,5.0);
  fHist_Deuteron_TPC_Chi2overNDF->GetXaxis()->SetTitle("#it{#chi}^{2}/ndf");
  fHist_Deuteron_TPC_Chi2overNDF->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_Chi2overNDF);

  fHist_Deuteron_ITS_nCluster = new TH1F("fHist_Deuteron_ITS_nCluster","Number of ITS clusters of deuteron tracks",20,0.0,20.0);
  fHist_Deuteron_ITS_nCluster->GetXaxis()->SetTitle("n_{Cluster} ITS");
  fHist_Deuteron_ITS_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_ITS_nCluster);

  fHist_Deuteron_TOF_Signal = new TH1F("fHist_Deuteron_TOF_Signal","TOF signal of deuteron tracks",200,0.0,200.0);
  fHist_Deuteron_TOF_Signal->GetXaxis()->SetTitle("TOF signal (ns)");
  fHist_Deuteron_TOF_Signal->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_Signal);

  fHist_Deuteron_DCAxy = new TH2F("fHist_Deuteron_DCAxy","DCA_{xy} distribution of deuteron tracks",240,0.0,6.0,1000,-5.0,+5.0);
  fHist_Deuteron_DCAxy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_DCAxy->GetYaxis()->SetTitle("DCA_{xy} (cm)");
  fHistList_Deuteron->Add(fHist_Deuteron_DCAxy);

  fHist_Deuteron_DCAz = new TH2F("fHist_Deuteron_DCAz","DCA_{z} distribution of deuteron tracks",240,0.0,6.0,1000,-5.0,+5.0);
  fHist_Deuteron_DCAz->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_DCAz->GetYaxis()->SetTitle("DCA_{z} (cm)");
  fHistList_Deuteron->Add(fHist_Deuteron_DCAz);

  fHist_Deuteron_TPC_nSigma_pT = new TH2F("fHist_Deuteron_TPC_nSigma_pT","n#it{#sigma}_{TPC} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_TPC_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_TPC_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_nSigma_pT);

  fHist_Deuteron_TPC_nSigma_p = new TH2F("fHist_Deuteron_TPC_nSigma_p","n#it{#sigma}_{TPC} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_TPC_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_TPC_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_nSigma_p);

  fHist_Deuteron_TOF_nSigma_pT = new TH2F("fHist_Deuteron_TOF_nSigma_pT","n#it{#sigma}_{TOF} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_TOF_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_TOF_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_nSigma_pT);

  fHist_Deuteron_TOF_nSigma_p = new TH2F("fHist_Deuteron_TOF_nSigma_p","n#it{#sigma}_{TOF} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_TOF_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_TOF_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_nSigma_p);

  fHist_Deuteron_TOF_nSigmaMassSq_pT = new TH2F("fHist_Deuteron_TOF_nSigmaMassSq_pT","n#it{#sigma}_{TOF} m^{2} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_TOF_nSigmaMassSq_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_TOF_nSigmaMassSq_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}(#it{m}^{2})");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_nSigmaMassSq_pT);

  fHist_Deuteron_TOF_nSigmaMassSq_p = new TH2F("fHist_Deuteron_TOF_nSigmaMassSq_p","n#it{#sigma}_{TOF} m^{2} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_TOF_nSigmaMassSq_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_TOF_nSigmaMassSq_p->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}(#it{m}^{2})");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_nSigmaMassSq_p);

  fHist_Deuteron_ITS_nSigma_pT = new TH2F("fHist_Deuteron_ITS_nSigma_pT","n#it{#sigma}_{ITS} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_ITS_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_ITS_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_Deuteron->Add(fHist_Deuteron_ITS_nSigma_pT);

  fHist_Deuteron_ITS_nSigma_p = new TH2F("fHist_Deuteron_ITS_nSigma_p","n#it{#sigma}_{ITS} of deuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_Deuteron_ITS_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_ITS_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_Deuteron->Add(fHist_Deuteron_ITS_nSigma_p);

  fHist_Deuteron_TPC_dEdx_pT = new TH2F("fHist_Deuteron_TPC_dEdx_pT","d#it{E}/d#it{x} of deuterons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_Deuteron_TPC_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_TPC_dEdx_pT->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_dEdx_pT);

  fHist_Deuteron_TPC_dEdx_p = new TH2F("fHist_Deuteron_TPC_dEdx_p","d#it{E}/d#it{x} of deuterons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_Deuteron_TPC_dEdx_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_TPC_dEdx_p->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_Deuteron->Add(fHist_Deuteron_TPC_dEdx_p);

  fHist_Deuteron_TOF_Beta_pT = new TH2F("fHist_Deuteron_TOF_Beta_pT","#it{#beta} of deuterons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_Deuteron_TOF_Beta_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_TOF_Beta_pT->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_Beta_pT);

  fHist_Deuteron_TOF_Beta_p = new TH2F("fHist_Deuteron_TOF_Beta_p","#it{#beta} of deuterons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_Deuteron_TOF_Beta_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_TOF_Beta_p->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_Beta_p);

  fHist_Deuteron_TOF_MassSquare_pT = new TH2F("fHist_Deuteron_TOF_MassSquare_pT","#it{m}^{2} of deuterons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_Deuteron_TOF_MassSquare_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_TOF_MassSquare_pT->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_MassSquare_pT);

  fHist_Deuteron_TOF_MassSquare_p = new TH2F("fHist_Deuteron_TOF_MassSquare_p","#it{m}^{2} of deuterons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_Deuteron_TOF_MassSquare_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_Deuteron_TOF_MassSquare_p->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_MassSquare_p);

  fHist_Deuteron_ITS_dEdx_pT = new TH2F("fHist_Deuteron_ITS_dEdx_pT","d#it{E}/d#it{x} of deuterons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_Deuteron_ITS_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_Deuteron_ITS_dEdx_pT->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_Deuteron->Add(fHist_Deuteron_ITS_dEdx_pT);

  fHist_Deuteron_ITS_dEdx_p = new TH2F("fHist_Deuteron_ITS_dEdx_p","d#it{E}/d#it{x} of deuterons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_Deuteron_ITS_dEdx_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_Deuteron_ITS_dEdx_p->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_Deuteron->Add(fHist_Deuteron_ITS_dEdx_p);




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for p-d pairs ++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_ProtonDeuteron_CutCounter = new TH1F("fHist_ProtonDeuteron_CutCounter","cut counter of p-d selection",3,0.0,3.0);
  fHist_ProtonDeuteron_CutCounter->GetYaxis()->SetTitle("counts");
  fHist_ProtonDeuteron_CutCounter->GetXaxis()->SetNdivisions(103,false);
  fHist_ProtonDeuteron_CutCounter->GetXaxis()->CenterLabels();
  fHist_ProtonDeuteron_CutCounter->GetXaxis()->SetBinLabel(1,"total pairs");
  fHist_ProtonDeuteron_CutCounter->GetXaxis()->SetBinLabel(2,"same ID");
  fHist_ProtonDeuteron_CutCounter->GetXaxis()->SetBinLabel(3,"close pair rejection");
  fHist_ProtonDeuteron_CutCounter->SetStats(0);
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_CutCounter);

  fHist_ProtonDeuteron_SED = new TH1F("fHist_ProtonDeuteron_SED","same-event distribution of p-d pairs",750,0.0,3.0);
  fHist_ProtonDeuteron_SED->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHist_ProtonDeuteron_SED->GetYaxis()->SetTitle("counts");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_SED);

  fHist_ProtonDeuteron_MED = new TH1F("fHist_ProtonDeuteron_MED","mixed-event distribution of p-d pairs",750,0.0,3.0);
  fHist_ProtonDeuteron_MED->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHist_ProtonDeuteron_MED->GetYaxis()->SetTitle("counts");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_MED);

  fHist_ProtonDeuteron_RPD = new TH1F("fHist_ProtonDeuteron_RPD","rotated pair distribution of p-d pairs",750,0.0,3.0);
  fHist_ProtonDeuteron_RPD->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHist_ProtonDeuteron_RPD->GetYaxis()->SetTitle("counts");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_RPD);

  fHist_ProtonDeuteron_PairsPerEvent = new TH1F("fHist_ProtonDeuteron_PairsPerEvent","number of selected p-d pairs per event",200,0.0,200);
  fHist_ProtonDeuteron_PairsPerEvent->GetXaxis()->SetTitle("#it{N}_{p-d} selected per event");
  fHist_ProtonDeuteron_PairsPerEvent->GetYaxis()->SetTitle("counts");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_PairsPerEvent);

  fHist_ProtonDeuteron_EventsForMixing = new TH1F("fHist_ProtonDeuteron_EventsForMixing","Event for p-d event-mixing found?",2,0.0,2.0);
  fHist_ProtonDeuteron_EventsForMixing->GetXaxis()->SetNdivisions(102,false);
  fHist_ProtonDeuteron_EventsForMixing->GetXaxis()->CenterLabels();
  fHist_ProtonDeuteron_EventsForMixing->GetXaxis()->SetBinLabel(1,"no event found");
  fHist_ProtonDeuteron_EventsForMixing->GetXaxis()->SetBinLabel(2,"events found");
  fHist_ProtonDeuteron_EventsForMixing->GetYaxis()->SetTitle("counts");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_EventsForMixing);

  fHist_ProtonDeuteron_AngleOfPairs = new TH2F("fHist_ProtonDeuteron_AngleOfPairs","angle between p-d pair",360,0.0,180,750,0.0,3.0);
  fHist_ProtonDeuteron_AngleOfPairs->GetXaxis()->SetTitle("angle (#circ)");
  fHist_ProtonDeuteron_AngleOfPairs->GetYaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_AngleOfPairs);

  fHist_ProtonDeuteron_PairMultiplicity = new TH2F("fHist_ProtonDeuteron_PairMultiplicity","number of protons and deuterons per event",50,0,50,50,0,50);
  fHist_ProtonDeuteron_PairMultiplicity->GetXaxis()->SetTitle("number of deuterons");
  fHist_ProtonDeuteron_PairMultiplicity->GetYaxis()->SetTitle("number of protons");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_PairMultiplicity);

  fHist_ProtonDeuteron_pT = new TH2F("fHist_ProtonDeuteron_pT","pair #it{p}_{T}",240,0.0,6.0,240,0.0,6.0);
  fHist_ProtonDeuteron_pT->GetXaxis()->SetTitle("proton #it{p}_{T} (GeV/#it{c})");
  fHist_ProtonDeuteron_pT->GetYaxis()->SetTitle("deuteron #it{p}_{T} (GeV/#it{c})");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_pT);

  fHist_ProtonDeuteron_Eta = new TH2F("fHist_ProtonDeuteron_Eta","pair #it{#eta}",200,-1.5,1.5,200,-1.5,1.5);
  fHist_ProtonDeuteron_Eta->GetXaxis()->SetTitle("proton #it{#eta}");
  fHist_ProtonDeuteron_Eta->GetYaxis()->SetTitle("deuteron #it{#eta}");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_Eta);

  fHist_ProtonDeuteron_Centrality = new TH2F("fHist_ProtonDeuteron_Centrality","p-d centrality",100,0.0,100.0,750,0.0,3.0);
  fHist_ProtonDeuteron_Centrality->GetXaxis()->SetTitle("centrality (%)");
  fHist_ProtonDeuteron_Centrality->GetYaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_Centrality);

  fHist_ProtonDeuteron_VertexZ = new TH2F("fHist_ProtonDeuteron_VertexZ","p-d primary vertex z-position",250,-12.0,+12.0,750,0.0,3.0);
  fHist_ProtonDeuteron_VertexZ->GetXaxis()->SetTitle("z-position of primary vertex (cm)");
  fHist_ProtonDeuteron_VertexZ->GetYaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_VertexZ);

  double Centrality_min = 0.0;
  double Centrality_max = 0.0;

  if(fCollisionSystem == 1) Centrality_min = CentralityBins_Central[0];
  if(fCollisionSystem == 1) Centrality_max = CentralityBins_Central[nCentralityBins];
  if(fCollisionSystem == 2) Centrality_min = CentralityBins_SemiCentral[0];
  if(fCollisionSystem == 2) Centrality_max = CentralityBins_SemiCentral[nCentralityBins];

  fHist_ProtonDeuteron_UsedEventsInPool = new TH2F("fHist_ProtonDeuteron_UsedEventsInPool","used events in p-d event-mixing",nCentralityBins,Centrality_min,Centrality_max,nZvtxBins,ZvtxBins[0],ZvtxBins[nZvtxBins]);
  fHist_ProtonDeuteron_UsedEventsInPool->GetXaxis()->SetTitle("centrality (%)");
  fHist_ProtonDeuteron_UsedEventsInPool->GetYaxis()->SetTitle("z-vertex position (cm)");
  fHistList_ProtonDeuteron->Add(fHist_ProtonDeuteron_UsedEventsInPool);

  
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for antiprotons ++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_AntiProton_CutCounter = new TH1F("fHist_AntiProton_CutCounter","cut counter of antiproton selection",20,0.0,20.0);
  fHist_AntiProton_CutCounter->GetYaxis()->SetTitle("counts");
  fHist_AntiProton_CutCounter->GetXaxis()->SetNdivisions(120,false);
  fHist_AntiProton_CutCounter->GetXaxis()->CenterLabels();
  fHist_AntiProton_CutCounter->GetXaxis()->SetLabelSize(0.025);
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(1,"no cuts");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(2,"TPC PID status");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(3,"n#sigma_{TPC}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(4,"DCA_{xy}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(5,"DCA_{z}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(6,"TOF PID status");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(7,"n#sigma_{TOF}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(8,"reject smaller #sigma");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(9,"FilterBit");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(10,"#it{p}_{T}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(11,"charge");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(12,"#eta");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(13,"n_{Cluster (TPC)}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(14,"n_{CrossedRows (TPC)}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(15,"n_{SharedCluster (TPC)}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(16,"#frac{n_{CrossedRows}}{n_{FindableCluster}}(TPC)");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(17,"low-#it{p}_{T} #pi");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(18,"ITS PID status");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(19,"n#sigma_{ITS}");
  fHist_AntiProton_CutCounter->GetXaxis()->SetBinLabel(20,"n_{Cluster} (ITS)");
  fHist_AntiProton_CutCounter->SetStats(0);
  fHistList_AntiProton->Add(fHist_AntiProton_CutCounter);

  fHist_AntiProton_pT = new TH1F("fHist_AntiProton_pT","#it{p}_{T} distribution of antiprotons",240,0.0,6.0);
  fHist_AntiProton_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_pT->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_pT);

  fHist_AntiProton_p = new TH1F("fHist_AntiProton_p","#it{p} distribution of antiprotons",240,0.0,6.0);
  fHist_AntiProton_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_AntiProton_p->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_p);

  fHist_AntiProton_pTPC = new TH1F("fHist_AntiProton_pTPC","#it{p}_{TPC} distribution of antiprotons",240,0.0,6.0);
  fHist_AntiProton_pTPC->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_pTPC->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_pTPC);

  fHist_AntiProton_Eta = new TH1F("fHist_AntiProton_Eta","#it{#eta} distribution of antiprotons",200,-1.5,+1.5);
  fHist_AntiProton_Eta->GetXaxis()->SetTitle("#it{#eta}");
  fHist_AntiProton_Eta->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_Eta);

  fHist_AntiProton_Phi = new TH1F("fHist_AntiProton_Phi","#it{#phi} distribution of antiprotons",360,0.0,360.0);
  fHist_AntiProton_Phi->GetXaxis()->SetTitle("#it{#phi} (#circ)");
  fHist_AntiProton_Phi->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_Phi);

  fHist_AntiProton_TPC_nCluster = new TH1F("fHist_AntiProton_TPC_nCluster","Number of TPC clusters of antiproton tracks",160,0.0,160.0);
  fHist_AntiProton_TPC_nCluster->GetXaxis()->SetTitle("n_{Cluster,TPC}");
  fHist_AntiProton_TPC_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_nCluster);

  fHist_AntiProton_TPC_CrossedRows = new TH1F("fHist_AntiProton_TPC_CrossedRows","Number of crossed TPC rows of antiproton tracks",160,0.0,160.0);
  fHist_AntiProton_TPC_CrossedRows->GetXaxis()->SetTitle("n_{Rows,TPC}");
  fHist_AntiProton_TPC_CrossedRows->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_CrossedRows);

  fHist_AntiProton_TPC_RatioRowsCluster = new TH1F("fHist_AntiProton_TPC_RatioRowsCluster","Crossed TPC rows / Findable clusters for antiprotons",150,0.0,1.5);
  fHist_AntiProton_TPC_RatioRowsCluster->GetXaxis()->SetTitle("n_{Rows,TPC}/n_{Cluster,TPC}");
  fHist_AntiProton_TPC_RatioRowsCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_RatioRowsCluster);

  fHist_AntiProton_TPC_SharedCluster = new TH1F("fHist_AntiProton_TPC_SharedCluster","Number of shared TPC clusters of antiproton tracks",160,0.0,160.0);
  fHist_AntiProton_TPC_SharedCluster->GetXaxis()->SetTitle("n_{Shared Clusters,TPC}");
  fHist_AntiProton_TPC_SharedCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_SharedCluster);

  fHist_AntiProton_TPC_Chi2overNDF = new TH1F("fHist_AntiProton_TPC_Chi2overNDF","#it{#chi}^{2}/ndf of TPC antiproton tracks",100,0.0,5.0);
  fHist_AntiProton_TPC_Chi2overNDF->GetXaxis()->SetTitle("#it{#chi}^{2}/ndf");
  fHist_AntiProton_TPC_Chi2overNDF->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_Chi2overNDF);

  fHist_AntiProton_ITS_nCluster = new TH1F("fHist_AntiProton_ITS_nCluster","Number of ITS clusters of antiproton tracks",20,0.0,20.0);
  fHist_AntiProton_ITS_nCluster->GetXaxis()->SetTitle("n_{Cluster} ITS");
  fHist_AntiProton_ITS_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_ITS_nCluster);

  fHist_AntiProton_TOF_Signal = new TH1F("fHist_AntiProton_TOF_Signal","TOF signal of antiproton tracks",200,0.0,200.0);
  fHist_AntiProton_TOF_Signal->GetXaxis()->SetTitle("TOF signal (ns)");
  fHist_AntiProton_TOF_Signal->GetYaxis()->SetTitle("counts");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_Signal);

  fHist_AntiProton_DCAxy = new TH2F("fHist_AntiProton_DCAxy","DCA_{xy} distribution of antiproton tracks",240,0.0,6.0,1000,-3.0,+3.0);
  fHist_AntiProton_DCAxy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_DCAxy->GetYaxis()->SetTitle("DCA_{xy} (cm)");
  fHistList_AntiProton->Add(fHist_AntiProton_DCAxy);

  fHist_AntiProton_DCAz = new TH2F("fHist_AntiProton_DCAz","DCA_{z} distribution of antiproton tracks",240,0.0,6.0,1000,-3.0,+3.0);
  fHist_AntiProton_DCAz->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_DCAz->GetYaxis()->SetTitle("DCA_{z} (cm)");
  fHistList_AntiProton->Add(fHist_AntiProton_DCAz);

  fHist_AntiProton_TPC_nSigma_pT = new TH2F("fHist_AntiProton_TPC_nSigma_pT","n#it{#sigma}_{TPC} of antiproton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiProton_TPC_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_TPC_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_nSigma_pT);

  fHist_AntiProton_TPC_nSigma_p = new TH2F("fHist_AntiProton_TPC_nSigma_p","n#it{#sigma}_{TPC} of antiproton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiProton_TPC_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_TPC_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_nSigma_p);

  fHist_AntiProton_TOF_nSigma_pT = new TH2F("fHist_AntiProton_TOF_nSigma_pT","n#it{#sigma}_{TOF} of antiproton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiProton_TOF_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_TOF_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_nSigma_pT);

  fHist_AntiProton_TOF_nSigma_p = new TH2F("fHist_AntiProton_TOF_nSigma_p","n#it{#sigma}_{TOF} of antiproton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiProton_TOF_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_TOF_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_nSigma_p);

  fHist_AntiProton_ITS_nSigma_pT = new TH2F("fHist_AntiProton_ITS_nSigma_pT","n#it{#sigma}_{ITS} of antiproton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiProton_ITS_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_ITS_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_AntiProton->Add(fHist_AntiProton_ITS_nSigma_pT);

  fHist_AntiProton_ITS_nSigma_p = new TH2F("fHist_AntiProton_ITS_nSigma_p","n#it{#sigma}_{ITS} of antiproton tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiProton_ITS_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_ITS_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_AntiProton->Add(fHist_AntiProton_ITS_nSigma_p);

  fHist_AntiProton_TPC_dEdx_pT = new TH2F("fHist_AntiProton_TPC_dEdx_pT","d#it{E}/d#it{x} of antiprotons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiProton_TPC_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_TPC_dEdx_pT->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_dEdx_pT);

  fHist_AntiProton_TPC_dEdx_p = new TH2F("fHist_AntiProton_TPC_dEdx_p","d#it{E}/d#it{x} of antiprotons meausred in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiProton_TPC_dEdx_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_TPC_dEdx_p->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_AntiProton->Add(fHist_AntiProton_TPC_dEdx_p);

  fHist_AntiProton_TOF_Beta_pT = new TH2F("fHist_AntiProton_TOF_Beta_pT","#it{#beta} of antiprotons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_AntiProton_TOF_Beta_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_TOF_Beta_pT->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_Beta_pT);

  fHist_AntiProton_TOF_Beta_p = new TH2F("fHist_AntiProton_TOF_Beta_p","#it{#beta} of antiprotons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_AntiProton_TOF_Beta_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_TOF_Beta_p->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_Beta_p);

  fHist_AntiProton_TOF_MassSquare_pT = new TH2F("fHist_AntiProton_TOF_MassSquare_pT","#it{m}^{2} of antiprotons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_AntiProton_TOF_MassSquare_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_TOF_MassSquare_pT->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_MassSquare_pT);

  fHist_AntiProton_TOF_MassSquare_p = new TH2F("fHist_AntiProton_TOF_MassSquare_p","#it{m}^{2} of antiprotons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_AntiProton_TOF_MassSquare_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiProton_TOF_MassSquare_p->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_AntiProton->Add(fHist_AntiProton_TOF_MassSquare_p);

  fHist_AntiProton_ITS_dEdx_pT = new TH2F("fHist_AntiProton_ITS_dEdx_pT","d#it{E}/d#it{x} of antiprotons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiProton_ITS_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiProton_ITS_dEdx_pT->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_AntiProton->Add(fHist_AntiProton_ITS_dEdx_pT);

  fHist_AntiProton_ITS_dEdx_p = new TH2F("fHist_AntiProton_ITS_dEdx_p","d#it{E}/d#it{x} of antiprotons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiProton_ITS_dEdx_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_AntiProton_ITS_dEdx_p->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_AntiProton->Add(fHist_AntiProton_ITS_dEdx_p);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for antideuterons ++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_AntiDeuteron_CutCounter = new TH1F("fHist_AntiDeuteron_CutCounter","cut counter of antideuteron selection",21,0.0,21.0);
  fHist_AntiDeuteron_CutCounter->GetYaxis()->SetTitle("counts");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetNdivisions(121,false);
  fHist_AntiDeuteron_CutCounter->GetXaxis()->CenterLabels();
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetLabelSize(0.025);
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(1,"no cuts");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(2,"TPC PID status");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(3,"n#sigma_{TPC}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(4,"DCA_{xy}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(5,"DCA_{z}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(6,"TOF PID status");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(7,"n#sigma_{TOF}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(8,"reject smaller #sigma");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(9,"n#sigma_{TOF} m^{2}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(10,"FilterBit");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(11,"#it{p}_{T}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(12,"charge");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(13,"#eta");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(14,"n_{Cluster (TPC)}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(15,"n_{CrossedRows (TPC)}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(16,"n_{SharedCluster (TPC)}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(17,"#frac{n_{CrossedRows}}{n_{FindableCluster}}(TPC)");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(18,"low-#it{p}_{T} #pi");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(19,"ITS PID status");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(20,"n#sigma_{ITS}");
  fHist_AntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(21,"n_{Cluster} (ITS)");
  fHist_AntiDeuteron_CutCounter->SetStats(0);
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_CutCounter);

  fHist_AntiDeuteron_pT = new TH1F("fHist_AntiDeuteron_pT","#it{p}_{T} distribution of antideuterons",240,0.0,6.0);
  fHist_AntiDeuteron_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_pT->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_pT);

  fHist_AntiDeuteron_p = new TH1F("fHist_AntiDeuteron_p","#it{p} distribution of antideuterons",240,0.0,6.0);
  fHist_AntiDeuteron_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_AntiDeuteron_p->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_p);

  fHist_AntiDeuteron_pTPC = new TH1F("fHist_AntiDeuteron_pTPC","#it{p}_{TPC} distribution of antideuterons",240,0.0,6.0);
  fHist_AntiDeuteron_pTPC->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_pTPC->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_pTPC);

  fHist_AntiDeuteron_Eta = new TH1F("fHist_AntiDeuteron_Eta","#it{#eta} distribution of antideuterons",200,-1.5,+1.5);
  fHist_AntiDeuteron_Eta->GetXaxis()->SetTitle("#it{#eta}");
  fHist_AntiDeuteron_Eta->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_Eta);

  fHist_AntiDeuteron_Phi = new TH1F("fHist_AntiDeuteron_Phi","#it{#phi} distribution of antideuterons",360,0.0,360.0);
  fHist_AntiDeuteron_Phi->GetXaxis()->SetTitle("#it{#phi} (#circ)");
  fHist_AntiDeuteron_Phi->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_Phi);

  fHist_AntiDeuteron_TPC_nCluster = new TH1F("fHist_AntiDeuteron_TPC_nCluster","Number of TPC clusters of antideuteron tracks",160,0.0,160.0);
  fHist_AntiDeuteron_TPC_nCluster->GetXaxis()->SetTitle("n_{Cluster,TPC}");
  fHist_AntiDeuteron_TPC_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_nCluster);

  fHist_AntiDeuteron_TPC_CrossedRows = new TH1F("fHist_AntiDeuteron_TPC_CrossedRows","Number of crossed TPC rows of antideuteron tracks",160,0.0,160.0);
  fHist_AntiDeuteron_TPC_CrossedRows->GetXaxis()->SetTitle("n_{Rows,TPC}");
  fHist_AntiDeuteron_TPC_CrossedRows->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_CrossedRows);

  fHist_AntiDeuteron_TPC_RatioRowsCluster = new TH1F("fHist_AntiDeuteron_TPC_RatioRowsCluster","Crossed TPC rows / Findable clusters for antideuterons",150,0.0,1.5);
  fHist_AntiDeuteron_TPC_RatioRowsCluster->GetXaxis()->SetTitle("n_{Rows,TPC}/n_{Cluster,TPC}");
  fHist_AntiDeuteron_TPC_RatioRowsCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_RatioRowsCluster);

  fHist_AntiDeuteron_TPC_SharedCluster = new TH1F("fHist_AntiDeuteron_TPC_SharedCluster","Number of shared TPC clusters of antideuteron tracks",160,0.0,160.0);
  fHist_AntiDeuteron_TPC_SharedCluster->GetXaxis()->SetTitle("n_{Shared Clusters,TPC}");
  fHist_AntiDeuteron_TPC_SharedCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_SharedCluster);

  fHist_AntiDeuteron_TPC_Chi2overNDF = new TH1F("fHist_AntiDeuteron_TPC_Chi2overNDF","#it{#chi}^{2}/ndf of TPC antideuteron tracks",100,0.0,5.0);
  fHist_AntiDeuteron_TPC_Chi2overNDF->GetXaxis()->SetTitle("#it{#chi}^{2}/ndf");
  fHist_AntiDeuteron_TPC_Chi2overNDF->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_Chi2overNDF);

  fHist_AntiDeuteron_ITS_nCluster = new TH1F("fHist_AntiDeuteron_ITS_nCluster","Number of ITS clusters of antideuteron tracks",20,0.0,20.0);
  fHist_AntiDeuteron_ITS_nCluster->GetXaxis()->SetTitle("n_{Cluster} ITS");
  fHist_AntiDeuteron_ITS_nCluster->GetYaxis()->SetTitle("counts");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_ITS_nCluster);

  fHist_AntiDeuteron_TOF_Signal = new TH1F("fHist_AntiDeuteron_TOF_Signal","TOF signal of antideuteron tracks",200,0.0,200.0);
  fHist_AntiDeuteron_TOF_Signal->GetXaxis()->SetTitle("TOF signal (ns)");
  fHist_AntiDeuteron_TOF_Signal->GetYaxis()->SetTitle("counts");
  fHistList_Deuteron->Add(fHist_Deuteron_TOF_Signal);

  fHist_AntiDeuteron_DCAxy = new TH2F("fHist_AntiDeuteron_DCAxy","DCA_{xy} distribution of antideuteron tracks",240,0.0,6.0,1000,-5.0,+5.0);
  fHist_AntiDeuteron_DCAxy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_DCAxy->GetYaxis()->SetTitle("DCA_{xy} (cm)");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_DCAxy);

  fHist_AntiDeuteron_DCAz = new TH2F("fHist_AntiDeuteron_DCAz","DCA_{z} distribution of antideuteron tracks",240,0.0,6.0,1000,-5.0,+5.0);
  fHist_AntiDeuteron_DCAz->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_DCAz->GetYaxis()->SetTitle("DCA_{z} (cm)");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_DCAz);

  fHist_AntiDeuteron_TPC_nSigma_pT = new TH2F("fHist_AntiDeuteron_TPC_nSigma_pT","n#it{#sigma}_{TPC} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_TPC_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_TPC_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_nSigma_pT);

  fHist_AntiDeuteron_TPC_nSigma_p = new TH2F("fHist_AntiDeuteron_TPC_nSigma_p","n#it{#sigma}_{TPC} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_TPC_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_TPC_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_nSigma_p);

  fHist_AntiDeuteron_TOF_nSigma_pT = new TH2F("fHist_AntiDeuteron_TOF_nSigma_pT","n#it{#sigma}_{TOF} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_TOF_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_nSigma_pT);

  fHist_AntiDeuteron_TOF_nSigma_p = new TH2F("fHist_AntiDeuteron_TOF_nSigma_p","n#it{#sigma}_{TOF} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_TOF_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_nSigma_p);

  fHist_AntiDeuteron_TOF_nSigmaMassSq_pT = new TH2F("fHist_AntiDeuteron_TOF_nSigmaMassSq_pT","n#it{#sigma}_{TOF} m^{2} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_TOF_nSigmaMassSq_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_nSigmaMassSq_pT->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}(#it{m}^{2})");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_nSigmaMassSq_pT);

  fHist_AntiDeuteron_TOF_nSigmaMassSq_p = new TH2F("fHist_AntiDeuteron_TOF_nSigmaMassSq_p","n#it{#sigma}_{TOF} m^{2} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_TOF_nSigmaMassSq_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_nSigmaMassSq_p->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}(#it{m}^{2})");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_nSigmaMassSq_p);

  fHist_AntiDeuteron_ITS_nSigma_pT = new TH2F("fHist_AntiDeuteron_ITS_nSigma_pT","n#it{#sigma}_{ITS} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_ITS_nSigma_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_ITS_nSigma_pT->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_ITS_nSigma_pT);

  fHist_AntiDeuteron_ITS_nSigma_p = new TH2F("fHist_AntiDeuteron_ITS_nSigma_p","n#it{#sigma}_{ITS} of antideuteron tracks",240,0.0,6.0,1200,-60.0,+60.0);
  fHist_AntiDeuteron_ITS_nSigma_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_ITS_nSigma_p->GetYaxis()->SetTitle("n#it{#sigma}_{ITS}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_ITS_nSigma_p);

  fHist_AntiDeuteron_TPC_dEdx_pT = new TH2F("fHist_AntiDeuteron_TPC_dEdx_pT","d#it{E}/d#it{x} of antideuterons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiDeuteron_TPC_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_TPC_dEdx_pT->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_dEdx_pT);

  fHist_AntiDeuteron_TPC_dEdx_p = new TH2F("fHist_AntiDeuteron_TPC_dEdx_p","d#it{E}/d#it{x} of antideuterons measured in TPC",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiDeuteron_TPC_dEdx_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_TPC_dEdx_p->GetYaxis()->SetTitle("TPC signal (a.u.)");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TPC_dEdx_p);

  fHist_AntiDeuteron_TOF_Beta_pT = new TH2F("fHist_AntiDeuteron_TOF_Beta_pT","#it{#beta} of antideuterons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_AntiDeuteron_TOF_Beta_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_Beta_pT->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_Beta_pT);

  fHist_AntiDeuteron_TOF_Beta_p = new TH2F("fHist_AntiDeuteron_TOF_Beta_p","#it{#beta} of antideuterons",240,0.0,6.0,1200,-0.1,+1.1);
  fHist_AntiDeuteron_TOF_Beta_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_Beta_p->GetYaxis()->SetTitle("#it{#beta}_{TOF}");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_Beta_p);

  fHist_AntiDeuteron_TOF_MassSquare_pT = new TH2F("fHist_AntiDeuteron_TOF_MassSquare_pT","#it{m}^{2} of antideuterons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_AntiDeuteron_TOF_MassSquare_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_MassSquare_pT->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_MassSquare_pT);

  fHist_AntiDeuteron_TOF_MassSquare_p = new TH2F("fHist_AntiDeuteron_TOF_MassSquare_p","#it{m}^{2} of antideuterons measured in TOF",240,0.0,6.0,1200,0.0,8.0);
  fHist_AntiDeuteron_TOF_MassSquare_p->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fHist_AntiDeuteron_TOF_MassSquare_p->GetYaxis()->SetTitle("#it{m}^{2}_{TOF} (GeV/#it{c}^{2})");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_TOF_MassSquare_p);

  fHist_AntiDeuteron_ITS_dEdx_pT = new TH2F("fHist_AntiDeuteron_ITS_dEdx_pT","d#it{E}/d#it{x} of antideuterons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiDeuteron_ITS_dEdx_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHist_AntiDeuteron_ITS_dEdx_pT->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_ITS_dEdx_pT);

  fHist_AntiDeuteron_ITS_dEdx_p = new TH2F("fHist_AntiDeuteron_ITS_dEdx_p","d#it{E}/d#it{x} of antideuterons measured in ITS",240,0.0,6.0,500,0.0,500.0);
  fHist_AntiDeuteron_ITS_dEdx_p->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fHist_AntiDeuteron_ITS_dEdx_p->GetYaxis()->SetTitle("ITS signal (a.u.)");
  fHistList_AntiDeuteron->Add(fHist_AntiDeuteron_ITS_dEdx_p);




  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ histograms for ap-ad pairs ++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fHist_AntiProtonAntiDeuteron_CutCounter = new TH1F("fHist_AntiProtonAntiDeuteron_CutCounter","cut counter of ap-ad selection",3,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_CutCounter->GetYaxis()->SetTitle("counts");
  fHist_AntiProtonAntiDeuteron_CutCounter->GetXaxis()->SetNdivisions(103,false);
  fHist_AntiProtonAntiDeuteron_CutCounter->GetXaxis()->CenterLabels();
  fHist_AntiProtonAntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(1,"total pairs");
  fHist_AntiProtonAntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(2,"same ID");
  fHist_AntiProtonAntiDeuteron_CutCounter->GetXaxis()->SetBinLabel(3,"close pair rejection");
  fHist_AntiProtonAntiDeuteron_CutCounter->SetStats(0);
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_CutCounter);

  fHist_AntiProtonAntiDeuteron_SED = new TH1F("fHist_AntiProtonAntiDeuteron_SED","same-event distribution of ap-ad pairs",750,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_SED->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHist_AntiProtonAntiDeuteron_SED->GetYaxis()->SetTitle("counts");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_SED);

  fHist_AntiProtonAntiDeuteron_MED = new TH1F("fHist_AntiProtonAntiDeuteron_MED","mixed-event distribution of ap-ad pairs",750,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_MED->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHist_AntiProtonAntiDeuteron_MED->GetYaxis()->SetTitle("counts");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_MED);

  fHist_AntiProtonAntiDeuteron_RPD = new TH1F("fHist_AntiProtonAntiDeuteron_RPD","rotated pair distribution of ap-ad pairs",750,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_RPD->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHist_AntiProtonAntiDeuteron_RPD->GetYaxis()->SetTitle("counts");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_RPD);

  fHist_AntiProtonAntiDeuteron_PairsPerEvent = new TH1F("fHist_AntiProtonAntiDeuteron_PairsPerEvent","number of selected ap-ad pairs per event",200,0.0,200);
  fHist_AntiProtonAntiDeuteron_PairsPerEvent->GetXaxis()->SetTitle("#it{N}_{#bar{p}-#bar{d}} selected per event");
  fHist_AntiProtonAntiDeuteron_PairsPerEvent->GetYaxis()->SetTitle("counts");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_PairsPerEvent);

  fHist_AntiProtonAntiDeuteron_EventsForMixing = new TH1F("fHist_AntiProtonAntiDeuteron_EventsForMixing","Event for #bar{p}-#bar{d} event-mixing found?",2,0.0,2.0);
  fHist_AntiProtonAntiDeuteron_EventsForMixing->GetXaxis()->SetNdivisions(102,false);
  fHist_AntiProtonAntiDeuteron_EventsForMixing->GetXaxis()->CenterLabels();
  fHist_AntiProtonAntiDeuteron_EventsForMixing->GetXaxis()->SetBinLabel(1,"no event found");
  fHist_AntiProtonAntiDeuteron_EventsForMixing->GetXaxis()->SetBinLabel(2,"events found");
  fHist_AntiProtonAntiDeuteron_EventsForMixing->GetYaxis()->SetTitle("counts");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_EventsForMixing);

  fHist_AntiProtonAntiDeuteron_AngleOfPairs = new TH2F("fHist_AntiProtonAntiDeuteron_AngleOfPairs","angle between #bar{p}-#bar{d} pair",360,0.0,180,750,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_AngleOfPairs->GetXaxis()->SetTitle("angle (#circ)");
  fHist_AntiProtonAntiDeuteron_AngleOfPairs->GetYaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_AngleOfPairs);

  fHist_AntiProtonAntiDeuteron_PairMultiplicity = new TH2F("fHist_AntiProtonAntiDeuteron_PairMultiplicity","number of antiprotons and antideuterons per event",50,0,50,50,0,50);
  fHist_AntiProtonAntiDeuteron_PairMultiplicity->GetXaxis()->SetTitle("number of antideuterons");
  fHist_AntiProtonAntiDeuteron_PairMultiplicity->GetYaxis()->SetTitle("number of antiprotons");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_PairMultiplicity);

  fHist_AntiProtonAntiDeuteron_pT = new TH2F("fHist_AntiProtonAntiDeuteron_pT","pair #it{p}_{T}",240,0.0,6.0,240,0.0,6.0);
  fHist_AntiProtonAntiDeuteron_pT->GetXaxis()->SetTitle("antiproton #it{p}_{T} (GeV/#it{c})");
  fHist_AntiProtonAntiDeuteron_pT->GetYaxis()->SetTitle("antideuteron #it{p}_{T} (GeV/#it{c})");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_pT);

  fHist_AntiProtonAntiDeuteron_Eta = new TH2F("fHist_AntiProtonAntiDeuteron_Eta","pair #it{#eta}",200,-1.5,1.5,200,-1.5,1.5);
  fHist_AntiProtonAntiDeuteron_Eta->GetXaxis()->SetTitle("antiproton #it{#eta}");
  fHist_AntiProtonAntiDeuteron_Eta->GetYaxis()->SetTitle("antideuteron #it{#eta}");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_Eta);

  fHist_AntiProtonAntiDeuteron_Centrality = new TH2F("fHist_AntiProtonAntiDeuteron_Centrality","#bar{p}-#bar{d} centrality",100,0.0,100.0,750,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_Centrality->GetXaxis()->SetTitle("centrality (%)");
  fHist_AntiProtonAntiDeuteron_Centrality->GetYaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_Centrality);

  fHist_AntiProtonAntiDeuteron_VertexZ = new TH2F("fHist_AntiProtonAntiDeuteron_VertexZ","#bar{p}-#bar{d} primary vertex z-position",250,-12.0,+12.0,750,0.0,3.0);
  fHist_AntiProtonAntiDeuteron_VertexZ->GetXaxis()->SetTitle("z-position of primary vertex (cm)");
  fHist_AntiProtonAntiDeuteron_VertexZ->GetYaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_VertexZ);

  fHist_AntiProtonAntiDeuteron_UsedEventsInPool = new TH2F("fHist_AntiProtonAntiDeuteron_UsedEventsInPool","used events in #bar{p}-#bar{d} event-mixing",nCentralityBins,Centrality_min,Centrality_max,nZvtxBins,ZvtxBins[0],ZvtxBins[nZvtxBins]);
  fHist_AntiProtonAntiDeuteron_UsedEventsInPool->GetXaxis()->SetTitle("centrality (%)");
  fHist_AntiProtonAntiDeuteron_UsedEventsInPool->GetYaxis()->SetTitle("z-vertex position (cm)");
  fHistList_AntiProtonAntiDeuteron->Add(fHist_AntiProtonAntiDeuteron_UsedEventsInPool);






/*

  // save some parameters of particles going into a low k*-pair
  fTree = new TTree("fTree","fTree");
  fTree->Branch("pTDeuteron",&DeuteronpT,"pTDeuteron/D");
  fTree->Branch("pTProton",&ProtonpT,"pTProton/D");
  fTree->Branch("pTLambda",&LambdapT,"pTLambda/D");
  fTree->Branch("pTPair",&PairpT,"pTPair/D");
  fTree->Branch("kstar",&kstar,"kstar/D");
  fTree->Branch("RunNumber",&RunNumber,"RunNumber/I")
  
*/

  // save particle labels
  ProtonTrackArray	  = new std::vector<int>;
  DeuteronTrackArray	  = new std::vector<int>;
  Lambdav0Array		  = new std::vector<int>;
  AntiProtonTrackArray	  = new std::vector<int>;
  AntiDeuteronTrackArray  = new std::vector<int>;
  AntiLambdav0Array	  = new std::vector<int>;





  PostData(1,fHistList_Event);
  PostData(2,fHistList_Proton);
  PostData(3,fHistList_Deuteron);
  PostData(4,fHistList_ProtonDeuteron);
  PostData(5,fHistList_AntiProton);
  PostData(6,fHistList_AntiDeuteron);
  PostData(7,fHistList_AntiProtonAntiDeuteron);



} // end of UserCreateOutputObjects









void AliAnalysisTask_pdLd::UserExec(Option_t*)
{

  bool debug = true;

  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  bool isMC = false;

//  AliMCEvent* fMCEvent = eventHandler->MCEvent();
//  if(!fMCEvent && isMC)::Fatal("AliAnalysisTask_pdLd::UserExec","No MC event found!");

  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAODEvent)::Fatal("AliAnalysisTask_pdLd::UserExec","No AOD event found!");

  fHeader = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  if(!fHeader)::Fatal("AliAnalysisTask_pdLd::UserExec","No Header found!");

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse)::Fatal("AliAnalysisTask_pdLd::UserExec","No PIDResponse found!");


  // debug analysis
  bool DebugAll		    = false;  // turn on all debugers
  bool DebugEventSelection  = false;
  bool DebugEventMixing	    = false;
  bool DebugPairSelection   = false;
  bool DebugPairRotation    = false;
  bool DebugRandomization   = false;

  if(DebugAll)
  {

    DebugEventMixing	= true;
    DebugEventSelection = true;
    DebugPairSelection	= true;
    DebugPairRotation	= true;

  } // end of DebugAll

 
//  bool isMC = false;
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
  if(!PrimaryVertex)::Warning("AliAnalsisTask_pdLd::UserExec","No AliAODVertex object found!");
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


  // fill event histograms
  fHist_Event_CutCounter->Fill(0);
  fHist_Event_Centrality->Fill(Centrality);
  fHist_Event_PrimVertexZ->Fill(PrimaryVertexZ);

  // fill event information into pool
  AliEventPool *EventPool_Deuteron	= fEventPoolManager->GetEventPool(Centrality,PrimaryVertexZ);
  AliEventPool *EventPool_AntiDeuteron	= fEventPoolManager->GetEventPool(Centrality,PrimaryVertexZ);

  // get event information
  int PeriodNumber	= fAODEvent->GetPeriodNumber();
  int RunNumber		= fAODEvent->GetRunNumber();
  int OrbitNumber	= fAODEvent->GetOrbitNumber();
  int BunchCrossNumber	= fAODEvent->GetBunchCrossNumber();
  int TimeStamp		= fAODEvent->GetTimeStamp();
  double bfield		= fAODEvent->GetMagneticField();
  

  // print event information
  if(DebugEventSelection)
  {

    cout << "" << endl;
    cout << "PeriodNumber:\t\t\t" << PeriodNumber << endl;
    cout << "RunNumber:\t\t\t" << RunNumber << endl;
    cout << "OrbitNumber:\t\t\t" << OrbitNumber << endl;
    cout << "BunchCrossNumber:\t\t" << BunchCrossNumber << endl;
    cout << "isPbPb:\t\t\t\t" << isPbPb << endl;
    cout << "Centrality:\t\t\t" << Centrality << " %" << endl;
    cout << "z-position of primary vertex:\t" << PrimaryVertexZ << " cm" << endl;
    cout << "Number of tracks in event:\t" << nTracks << endl;
    cout << "Magnetic field:\t\t\t" << bfield*0.1 << " T" << endl;
    cout << "Time stamp: " << TimeStamp << endl;

  }// end of DebugEventSelection





  // empty arrays containing labels of the previous event
  ProtonTrackArray->clear();
  DeuteronTrackArray->clear();
  Lambdav0Array->clear();
  AntiProtonTrackArray->clear();
  AntiDeuteronTrackArray->clear();
  AntiLambdav0Array->clear();


  // set up particle masses (GeV/c2)
  const double ProtonMass   = 0.9382720;
  const double DeuteronMass = 1.8756129;


  bool ProtonIsSelected	      = false;
  bool DeuteronIsSelected     = false;
  bool AntiProtonIsSelected   = false;
  bool AntiDeuteronIsSelected = false;

  int nRotatedPairs = 10;
  double kstar_max = 0.2;


  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ proton selection loop +++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// proton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pdLd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  
    fHist_Proton_CutCounter->Fill(0);

    // apply proton cuts
    bool PassedProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,true);
    if(!PassedProtonCuts) continue;

    double pT				= Track->Pt();
    double p				= Track->P();
    double pTPC				= Track->GetTPCmomentum();
    int nSharedClusterTPC		= Track->GetTPCnclsS();
    int nClusterTPCfindable		= Track->GetTPCNclsF();
    double nCrossedRowsTPC		= Track->GetTPCCrossedRows();
    double RatioRowsFindableClusterTPC	= ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);
      
    double nSigmaTOF  = -999.0;
    double beta	      = -999.0;
    double massSq     = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);
    if(statusTOF == AliPIDResponse::kDetPidOk)
    {

      nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
      beta	= CalculateBetaTOF(*Track);
      massSq	= CalculateMassSquareTOF(*Track);

    }

      
    double nSigmaITS = -999.0;
    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    if(statusITS == AliPIDResponse::kDetPidOk) nSigmaITS = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kProton);
 
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    // fill proton histograms

    // general
    fHist_Proton_pT->Fill(pT);
    fHist_Proton_p->Fill(p);
    fHist_Proton_pTPC->Fill(Track->GetTPCmomentum());
    fHist_Proton_Eta->Fill(Track->Eta());
    fHist_Proton_Phi->Fill(Track->Phi()*TMath::RadToDeg());
    fHist_Proton_DCAxy->Fill(pT,DCAxy);
    fHist_Proton_DCAz->Fill(pT,DCAz);

    // TPC
    fHist_Proton_TPC_CrossedRows->Fill(Track->GetTPCCrossedRows());
    fHist_Proton_TPC_RatioRowsCluster->Fill(RatioRowsFindableClusterTPC);
    fHist_Proton_TPC_SharedCluster->Fill(nSharedClusterTPC);
    fHist_Proton_TPC_Chi2overNDF->Fill(Track->GetTPCchi2perNDF());
    fHist_Proton_TPC_nCluster->Fill(Track->GetTPCNcls());
    fHist_Proton_TPC_nSigma_pT->Fill(pT,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton));
    fHist_Proton_TPC_nSigma_p->Fill(pTPC,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton));
    fHist_Proton_TPC_dEdx_pT->Fill(pT,Track->GetTPCsignal());
    fHist_Proton_TPC_dEdx_p->Fill(pTPC,Track->GetTPCsignal());

    // TOF
    fHist_Proton_TOF_Signal->Fill(Track->GetTOFsignal()/1000);
    fHist_Proton_TOF_nSigma_pT->Fill(pT,nSigmaTOF);
    fHist_Proton_TOF_nSigma_p->Fill(pTPC,nSigmaTOF);
    fHist_Proton_TOF_Beta_pT->Fill(pT,beta);
    fHist_Proton_TOF_Beta_p->Fill(pTPC,beta);
    fHist_Proton_TOF_MassSquare_pT->Fill(pT,massSq);
    fHist_Proton_TOF_MassSquare_p->Fill(pTPC,massSq);

    // ITS
    fHist_Proton_ITS_nCluster->Fill(Track->GetITSNcls());
    fHist_Proton_ITS_nSigma_pT->Fill(pT,nSigmaITS);
    fHist_Proton_ITS_nSigma_p->Fill(p,nSigmaITS);
    fHist_Proton_ITS_dEdx_pT->Fill(pT,Track->GetITSsignal());
    fHist_Proton_ITS_dEdx_p->Fill(p,Track->GetITSsignal());

 
    ProtonIsSelected = true;


    // save details of selected proton
    ProtonTrackArray->push_back(track);

  } // end of proton loop



  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ deuteron selection loop +++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// deuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pdLd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  
    fHist_Deuteron_CutCounter->Fill(0);

    // apply deuteron cuts
    bool PassedDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,true);
    if(!PassedDeuteronCuts) continue;

    double pT				= Track->Pt();
    double pTPC				= Track->GetTPCmomentum();
    int nSharedClusterTPC		= Track->GetTPCnclsS();
    int nClusterTPCfindable		= Track->GetTPCNclsF();
    double nCrossedRowsTPC = Track->GetTPCCrossedRows();
    double RatioRowsFindableClusterTPC	= ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);
    double p				= Track->P();

      
    double nSigmaTOF  = -999.0;
    double nSigmaTOFMassSq = -999.0; 
    double beta	      = -999.0;
    double massSq     = -999.0;


    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);

    if(statusTOF == AliPIDResponse::kDetPidOk)
    {

      nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
      beta	= CalculateBetaTOF(*Track);
      massSq	= CalculateMassSquareTOF(*Track);
      nSigmaTOFMassSq = CalculateDeuteronSigmaMassSquareTOF(pT,massSq,true);

    }

      
    double nSigmaITS = -999.0;
    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    if(statusITS == AliPIDResponse::kDetPidOk) nSigmaITS = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kDeuteron);

  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    // fill deuteron histograms

    // general
    fHist_Deuteron_pT->Fill(Track->Pt());
    fHist_Deuteron_p->Fill(Track->P());
    fHist_Deuteron_pTPC->Fill(pTPC);
    fHist_Deuteron_Eta->Fill(Track->Eta());
    fHist_Deuteron_Phi->Fill(Track->Phi()*TMath::RadToDeg());
    fHist_Deuteron_DCAxy->Fill(pT,DCAxy);
    fHist_Deuteron_DCAz->Fill(pT,DCAz);

    // TPC
    fHist_Deuteron_TPC_CrossedRows->Fill(Track->GetTPCCrossedRows());
    fHist_Deuteron_TPC_RatioRowsCluster->Fill(RatioRowsFindableClusterTPC);
    fHist_Deuteron_TPC_SharedCluster->Fill(nSharedClusterTPC);
    fHist_Deuteron_TPC_Chi2overNDF->Fill(Track->GetTPCchi2perNDF());
    fHist_Deuteron_TPC_nCluster->Fill(Track->GetTPCNcls());
    fHist_Deuteron_TPC_nSigma_pT->Fill(pT,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron));
    fHist_Deuteron_TPC_nSigma_p->Fill(pTPC,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron));
    fHist_Deuteron_TPC_dEdx_pT->Fill(pT,Track->GetTPCsignal());
    fHist_Deuteron_TPC_dEdx_p->Fill(pTPC,Track->GetTPCsignal());

    // TOF
    fHist_Deuteron_TOF_Signal->Fill(Track->GetTOFsignal()/1000);
    fHist_Deuteron_TOF_nSigma_pT->Fill(pT,nSigmaTOF);
    fHist_Deuteron_TOF_nSigma_p->Fill(pTPC,nSigmaTOF);
    fHist_Deuteron_TOF_nSigmaMassSq_pT->Fill(pT,nSigmaTOFMassSq);
    fHist_Deuteron_TOF_nSigmaMassSq_p->Fill(pTPC,nSigmaTOFMassSq);
    fHist_Deuteron_TOF_Beta_pT->Fill(pT,beta);
    fHist_Deuteron_TOF_Beta_p->Fill(pTPC,beta);
    fHist_Deuteron_TOF_MassSquare_pT->Fill(pT,massSq);
    fHist_Deuteron_TOF_MassSquare_p->Fill(pTPC,massSq);

    // ITS
    fHist_Deuteron_ITS_nCluster->Fill(Track->GetITSNcls());
    fHist_Deuteron_ITS_nSigma_pT->Fill(pT,nSigmaITS);
    fHist_Deuteron_ITS_nSigma_p->Fill(p,nSigmaITS);
    fHist_Deuteron_ITS_dEdx_pT->Fill(pT,Track->GetITSsignal());
    fHist_Deuteron_ITS_dEdx_p->Fill(p,Track->GetITSsignal());


    DeuteronIsSelected = true;

    // save details of selected deuteron
    DeuteronTrackArray->push_back(track);

  } // end of deuteron loop




  // number of events containing selected ...
  if(ProtonIsSelected)			      fHist_Event_CutCounter->Fill(1);
  if(DeuteronIsSelected)		      fHist_Event_CutCounter->Fill(2);
  if(ProtonIsSelected && DeuteronIsSelected)  fHist_Event_CutCounter->Fill(3);


  int nProtonsSelected = ProtonTrackArray->size();
  int nDeuteronsSelected = DeuteronTrackArray->size();
  int nProtonDeuteronPairs = 0;



  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ proton-deuteron pairing loop ++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(ProtonIsSelected && DeuteronIsSelected)
  {

    for(int Track1 = 0; Track1 < nProtonsSelected; Track1++)	// particle pair loop (proton loop)
    { 

      AliAODTrack *ProtonTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(ProtonTrackArray->at(Track1)));

      double MomentumProton[3];
      MomentumProton[0] = ProtonTrack->Px();
      MomentumProton[1] = ProtonTrack->Py();
      MomentumProton[2] = ProtonTrack->Pz();

      TLorentzVector LorentzVectorProton;
      LorentzVectorProton.SetXYZM(MomentumProton[0],MomentumProton[1],MomentumProton[2],ProtonMass);

      for(int Track2 = 0; Track2 < nDeuteronsSelected; Track2++) // particle pair loop (deuteron loop)
      {

	AliAODTrack *DeuteronTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(DeuteronTrackArray->at(Track2)));

	double MomentumDeuteron[3];
	MomentumDeuteron[0] = DeuteronTrack->Px();
	MomentumDeuteron[1] = DeuteronTrack->Py();
	MomentumDeuteron[2] = DeuteronTrack->Pz();

	TLorentzVector LorentzVectorDeuteron;
	LorentzVectorDeuteron.SetXYZM(MomentumDeuteron[0],MomentumDeuteron[1],MomentumDeuteron[2],DeuteronMass);

	fHist_ProtonDeuteron_CutCounter->Fill(0);

	int ProtonID = ProtonTrack->GetID();
	int DeuteronID = DeuteronTrack->GetID();

	if((ProtonID == DeuteronID) && DebugPairSelection)
	{
      
	  cout << "x-x-x-x-x-x-x-x-> same ID for proton and deuteron: " << ProtonID << endl;
	  cout << "Proton:\t\t\tDeuteron:" << endl;
	  cout << MomentumProton[0] << "\t\t\t" << MomentumDeuteron[0] << endl;
	  cout << MomentumProton[1] << "\t\t\t" << MomentumDeuteron[1] << endl;
	  cout << MomentumProton[2] << "\t\t\t" << MomentumDeuteron[2] << endl;

	}


	if(ProtonID == DeuteronID) continue;
	fHist_ProtonDeuteron_CutCounter->Fill(1);
	

	// apply eta and phi cut
	bool RejectClosePairs = ClosePairRejection(bfield,*ProtonTrack,*DeuteronTrack);
	if(RejectClosePairs) continue;
	fHist_ProtonDeuteron_CutCounter->Fill(2);

	double angle = LorentzVectorProton.Angle(LorentzVectorDeuteron.Vect());

	// do pair cleaning (avoid auto-correlations)
	// apply pair cuts

	TLorentzVector LorentzVectorPair = LorentzVectorProton + LorentzVectorDeuteron;
	double RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorProton,LorentzVectorDeuteron);

	if(DebugPairSelection) cout << "angle: " << angle*TMath::RadToDeg() << "°\t kstar: " << RelativeMomentum << endl;
  
	nProtonDeuteronPairs++;
	fHist_ProtonDeuteron_SED->Fill(RelativeMomentum);
	fHist_ProtonDeuteron_AngleOfPairs->Fill(angle*TMath::RadToDeg(),RelativeMomentum);
	fHist_ProtonDeuteron_Centrality->Fill(Centrality,RelativeMomentum);
	fHist_ProtonDeuteron_VertexZ->Fill(PrimaryVertexZ,RelativeMomentum);

	if(RelativeMomentum <= kstar_max){

	  fHist_ProtonDeuteron_pT->Fill(ProtonTrack->Pt(),DeuteronTrack->Pt());
	  fHist_ProtonDeuteron_Eta->Fill(ProtonTrack->Eta(),DeuteronTrack->Eta());

	}


	double phi = ProtonTrack->Phi();

	if(DebugPairRotation){

	  cout << "---------------------------------" << endl;
	  cout << "original phi of proton: " << phi << " rad" << endl;
	  cout << "original k* of p-d pair: " << RelativeMomentum << endl;
	  cout << "original Px of proton: " << MomentumProton[0] << endl;
	  cout << "original Py of proton: " << MomentumProton[1] << endl;
	  cout << "original Pz of proton: " << MomentumProton[2] << endl;
	  cout << "" << endl;

	} // end of DebugPairRotation

	// set up random generator
	TRandom3 *RandomGenerator = new TRandom3();
	RandomGenerator->SetSeed(0);

	// rotate pairs
	for(int iPair = 0; iPair < nRotatedPairs; iPair++)
	{

	  double RandomPhi = RandomGenerator->Uniform(0.0,2*TMath::Pi());
	  double phi_new = phi + RandomPhi;
	  ProtonTrack->SetPhi(phi_new);
      
	  double MomentumProtonRotated[3];
	  MomentumProtonRotated[0] = ProtonTrack->Px();
	  MomentumProtonRotated[1] = ProtonTrack->Py();
	  MomentumProtonRotated[2] = ProtonTrack->Pz();

	  TLorentzVector LorentzVectorProtonRotated; 
	  LorentzVectorProtonRotated.SetXYZM(MomentumProtonRotated[0],MomentumProtonRotated[1],MomentumProtonRotated[2],ProtonMass);

	  TLorentzVector LorentzVectorPairRotated = LorentzVectorProtonRotated + LorentzVectorDeuteron;
	  double RelativeMomentumRotated = CalculateRelativeMomentum(LorentzVectorPairRotated,LorentzVectorProtonRotated,LorentzVectorDeuteron);

	  fHist_ProtonDeuteron_RPD->Fill(RelativeMomentumRotated);

	  if(DebugPairRotation){

	    cout << "iPair: " << iPair << endl;
	    cout << "RandomPhi: " << RandomPhi << endl;
	    cout << "Phi_new: " << phi_new << endl;
	    cout << "Rotated k*: " << RelativeMomentumRotated << endl;
	    cout << "rotated Px of proton: " << MomentumProtonRotated[0] << endl;
	    cout << "rotated Py of proton: " << MomentumProtonRotated[1] << endl;
	    cout << "rotated Pz of proton: " << MomentumProtonRotated[2] << endl;

	  } // end of DebugPairRotation

	} // end of pair rotation

      } // end of particle pair loop (deuteron loop)

    } // end of particle pair loop (proton loop)

  } // end of if statement: ProtonIsSelected && DeuteronIsSelected


  if(DebugPairSelection)
  {

    if(!DebugEventSelection) cout << "" << endl;
    cout << "nProtons selected:\t\t" << nProtonsSelected << endl;
    cout << "nDeuterons selected:\t\t" << nDeuteronsSelected << endl;
    cout << "nPairs selected:\t\t" << nProtonDeuteronPairs << endl;

  } // end of DebugPairSelection

  fHist_ProtonDeuteron_PairsPerEvent->Fill(nProtonDeuteronPairs);


  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ proton-deuteron event mixing ++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(ProtonIsSelected && DeuteronIsSelected) // perform event-mixing
  {

    int nMixedEvents  = EventPool_Deuteron->GetCurrentNEvents();
    double ZvtxMin    = EventPool_Deuteron->GetZvtxMin();
    double ZvtxMax    = EventPool_Deuteron->GetZvtxMax();
    double MultMin    = EventPool_Deuteron->GetMultMin();
    double MultMax    = EventPool_Deuteron->GetMultMax();
    double PtMin      = EventPool_Deuteron->GetPtMin();
    double PtMax      = EventPool_Deuteron->GetPtMax();

    if(DebugEventMixing)
    {
    
      cout << "Number of events in buffer: " << nMixedEvents << endl;
      cout << "ZvtxMin: " << ZvtxMin << "\t\t ZvtxMax: " << ZvtxMax << endl;
      cout << "MultMin: " << MultMin << "\t\t MultMax: " << MultMax << endl;
      cout << "PtMin: " << PtMin << "\t\t PtMax: " << PtMax << endl;

    } // end of DebugEventMixing

    if(nMixedEvents == 0) fHist_ProtonDeuteron_EventsForMixing->Fill(0);
    if(nMixedEvents > 0)  fHist_ProtonDeuteron_EventsForMixing->Fill(1);

    for(int iMixedEvents = 0; iMixedEvents < nMixedEvents; iMixedEvents++){ // loop over events of the pool

      TObjArray *MixedEventArray = (TObjArray*)EventPool_Deuteron->GetEvent(iMixedEvents);
      if(!MixedEventArray) continue;

      AliAODTrackTiny *DeuteronTrack = (AliAODTrackTiny*) MixedEventArray->At(0);
      if(!DeuteronTrack) continue;

      fHist_ProtonDeuteron_UsedEventsInPool->Fill(DeuteronTrack->GetCentrality(),DeuteronTrack->GetPrimaryVertexZ());

      double MomentumDeuteron[3];
      MomentumDeuteron[0] = DeuteronTrack->Px();
      MomentumDeuteron[1] = DeuteronTrack->Py();
      MomentumDeuteron[2] = DeuteronTrack->Pz();
   
      TLorentzVector LorentzVectorDeuteron;
      LorentzVectorDeuteron.SetXYZM(MomentumDeuteron[0],MomentumDeuteron[1],MomentumDeuteron[2],DeuteronMass);

      for(int Track2 = 0; Track2 < nProtonsSelected; Track2++) // particle pair loop (proton loop)
      {

	AliAODTrack *ProtonTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(ProtonTrackArray->at(Track2)));
	if(!ProtonTrack) continue;

	double MomentumProton[3];
	MomentumProton[0] = ProtonTrack->Px();
	MomentumProton[1] = ProtonTrack->Py();
	MomentumProton[2] = ProtonTrack->Pz();

	TLorentzVector LorentzVectorProton;
	LorentzVectorProton.SetXYZM(MomentumProton[0],MomentumProton[1],MomentumProton[2],ProtonMass);

	TLorentzVector LorentzVectorPair = LorentzVectorProton + LorentzVectorDeuteron;
	double RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorProton,LorentzVectorDeuteron);

	if(TMath::IsNaN(RelativeMomentum)) continue;

	fHist_ProtonDeuteron_MED->Fill(RelativeMomentum);

      } // loop over protons

    } // loop over mixed events


    TObjArray *DeuteronObjectArray = new TObjArray();
    DeuteronObjectArray->SetOwner(true);

    int PickDeuteronNumber = 0;

    if(nDeuteronsSelected > 1)	// randomize deuteron selection if more than one deuteron is present in the current event
    {

      TRandom3 *RandomGenerator = new TRandom3();
      RandomGenerator->SetSeed(0);
      PickDeuteronNumber = RandomGenerator->Integer(nDeuteronsSelected); // produces a random integer between 0 and nDeuteronsSelected-1

    }

    if(DebugRandomization) std::cout << "nDeuteronsSelected: " << nDeuteronsSelected << "\t PickDeuteronNumber: " << PickDeuteronNumber << std::endl;

    AliAODTrack *DeuteronTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(DeuteronTrackArray->at(PickDeuteronNumber))); // array starts at 0
    if(!DeuteronTrack) std::cout << "AliAnalysisTask_pdLd: No DeuteronTrack received from the EventPool!" << std::endl;

    AliAODTrackTiny *TinyDeuteronTrack = new AliAODTrackTiny();
    TinyDeuteronTrack->InitFromTrack(DeuteronTrack,Centrality,PrimaryVertexZ);

    DeuteronObjectArray->Add(TinyDeuteronTrack);

    EventPool_Deuteron->UpdatePool(DeuteronObjectArray);

  } // end of proton and deuteron is selected








  int nPoolEvents = EventPool_Deuteron->GetCurrentNEvents();


  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antiproton selection loop +++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// antiproton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pdLd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  
    fHist_AntiProton_CutCounter->Fill(0);

    // apply antiproton cuts
    bool PassedAntiProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,false);
    if(!PassedAntiProtonCuts) continue;

    double pT				= Track->Pt();
    double p				= Track->P();
    double pTPC				= Track->GetTPCmomentum();
    int nSharedClusterTPC		= Track->GetTPCnclsS();
    int nClusterTPCfindable		= Track->GetTPCNclsF();
    double nCrossedRowsTPC		= Track->GetTPCCrossedRows();
    double RatioRowsFindableClusterTPC	= ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);

      
    double nSigmaTOF  = -999.0;
    double beta	      = -999.0;
    double massSq     = -999.0;

    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);

    if(statusTOF == AliPIDResponse::kDetPidOk)
    {

      nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);
      beta	= CalculateBetaTOF(*Track);
      massSq	= CalculateMassSquareTOF(*Track);

    }

      
    double nSigmaITS = -999.0;
    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    if(statusITS == AliPIDResponse::kDetPidOk) nSigmaITS = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kProton);
 
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    // fill antiproton histograms

    // general
    fHist_AntiProton_pT->Fill(pT);
    fHist_AntiProton_p->Fill(p);
    fHist_AntiProton_pTPC->Fill(Track->GetTPCmomentum());
    fHist_AntiProton_Eta->Fill(Track->Eta());
    fHist_AntiProton_Phi->Fill(Track->Phi()*TMath::RadToDeg());
    fHist_AntiProton_DCAxy->Fill(pT,DCAxy);
    fHist_AntiProton_DCAz->Fill(pT,DCAz);

    // TPC
    fHist_AntiProton_TPC_CrossedRows->Fill(Track->GetTPCCrossedRows());
    fHist_AntiProton_TPC_RatioRowsCluster->Fill(RatioRowsFindableClusterTPC);
    fHist_AntiProton_TPC_SharedCluster->Fill(nSharedClusterTPC);
    fHist_AntiProton_TPC_Chi2overNDF->Fill(Track->GetTPCchi2perNDF());
    fHist_AntiProton_TPC_nCluster->Fill(Track->GetTPCNcls());
    fHist_AntiProton_TPC_nSigma_pT->Fill(pT,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton));
    fHist_AntiProton_TPC_nSigma_p->Fill(pTPC,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton));
    fHist_AntiProton_TPC_dEdx_pT->Fill(pT,Track->GetTPCsignal());
    fHist_AntiProton_TPC_dEdx_p->Fill(pTPC,Track->GetTPCsignal());

    // TOF
    fHist_AntiProton_TOF_Signal->Fill(Track->GetTOFsignal()/1000);
    fHist_AntiProton_TOF_nSigma_pT->Fill(pT,nSigmaTOF);
    fHist_AntiProton_TOF_nSigma_p->Fill(pTPC,nSigmaTOF);
    fHist_AntiProton_TOF_Beta_pT->Fill(pT,beta);
    fHist_AntiProton_TOF_Beta_p->Fill(pTPC,beta);
    fHist_AntiProton_TOF_MassSquare_pT->Fill(pT,massSq);
    fHist_AntiProton_TOF_MassSquare_p->Fill(pTPC,massSq);

    // ITS
    fHist_AntiProton_ITS_nCluster->Fill(Track->GetITSNcls());
    fHist_AntiProton_ITS_nSigma_pT->Fill(pT,nSigmaITS);
    fHist_AntiProton_ITS_nSigma_p->Fill(p,nSigmaITS);
    fHist_AntiProton_ITS_dEdx_pT->Fill(pT,Track->GetITSsignal());
    fHist_AntiProton_ITS_dEdx_p->Fill(p,Track->GetITSsignal());

 
    AntiProtonIsSelected = true;

    // save details of selected proton
    AntiProtonTrackArray->push_back(track);

  } // end of antiproton loop








  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antideuteron selection loop +++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int track = 0; track < nTracks; track++)	// antideuteron loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pdLd::UserExec","No AliAODTrack found");
    if(!Track) continue;
  
    fHist_AntiDeuteron_CutCounter->Fill(0);

    // apply antideuteron cuts
    bool PassedAntiDeuteronCuts = CheckDeuteronCuts(*Track,*fPIDResponse,false);
    if(!PassedAntiDeuteronCuts) continue;

    double pT				= Track->Pt();
    double p				= Track->P();
    double pTPC				= Track->GetTPCmomentum();
    int nSharedClusterTPC		= Track->GetTPCnclsS();
    int nClusterTPCfindable		= Track->GetTPCNclsF();
    double nCrossedRowsTPC		= Track->GetTPCCrossedRows();
    double RatioRowsFindableClusterTPC	= ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);

      
    double nSigmaTOF	    = -999.0;
    double nSigmaTOFMassSq  = -999.0; 
    double beta		    = -999.0;
    double massSq	    = -999.0;


    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);

    if(statusTOF == AliPIDResponse::kDetPidOk)
    {

      nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kDeuteron);
      beta	= CalculateBetaTOF(*Track);
      massSq	= CalculateMassSquareTOF(*Track);
      nSigmaTOFMassSq = CalculateDeuteronSigmaMassSquareTOF(pT,massSq,false);

    }
      
    double nSigmaITS = -999.0;
    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    if(statusITS == AliPIDResponse::kDetPidOk) nSigmaITS = fPIDResponse->NumberOfSigmasITS(Track,AliPID::kDeuteron);
  
    float xv[2];
    float yv[3];
    Track->GetImpactParameters(xv,yv);
    double DCAxy = xv[0];
    double DCAz = xv[1];

    // fill antideuteron histograms

    // general
    fHist_AntiDeuteron_pT->Fill(Track->Pt());
    fHist_AntiDeuteron_p->Fill(Track->P());
    fHist_AntiDeuteron_pTPC->Fill(Track->GetTPCmomentum());
    fHist_AntiDeuteron_Eta->Fill(Track->Eta());
    fHist_AntiDeuteron_Phi->Fill(Track->Phi()*TMath::RadToDeg());
    fHist_AntiDeuteron_DCAxy->Fill(pT,DCAxy);
    fHist_AntiDeuteron_DCAz->Fill(pT,DCAz);

    // TPC
    fHist_AntiDeuteron_TPC_CrossedRows->Fill(Track->GetTPCCrossedRows());
    fHist_AntiDeuteron_TPC_RatioRowsCluster->Fill(RatioRowsFindableClusterTPC);
    fHist_AntiDeuteron_TPC_SharedCluster->Fill(nSharedClusterTPC);
    fHist_AntiDeuteron_TPC_Chi2overNDF->Fill(Track->GetTPCchi2perNDF());
    fHist_AntiDeuteron_TPC_nCluster->Fill(Track->GetTPCNcls());
    fHist_AntiDeuteron_TPC_nSigma_pT->Fill(pT,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron));
    fHist_AntiDeuteron_TPC_nSigma_p->Fill(pTPC,fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kDeuteron));
    fHist_AntiDeuteron_TPC_dEdx_pT->Fill(pT,Track->GetTPCsignal());
    fHist_AntiDeuteron_TPC_dEdx_p->Fill(pTPC,Track->GetTPCsignal());

    // TOF
    fHist_AntiDeuteron_TOF_Signal->Fill(Track->GetTOFsignal()/1000);
    fHist_AntiDeuteron_TOF_nSigma_pT->Fill(pT,nSigmaTOF);
    fHist_AntiDeuteron_TOF_nSigma_p->Fill(pTPC,nSigmaTOF);
    fHist_AntiDeuteron_TOF_nSigmaMassSq_pT->Fill(pT,nSigmaTOFMassSq);
    fHist_AntiDeuteron_TOF_nSigmaMassSq_p->Fill(pTPC,nSigmaTOFMassSq);
    fHist_AntiDeuteron_TOF_Beta_pT->Fill(pT,beta);
    fHist_AntiDeuteron_TOF_Beta_p->Fill(pTPC,beta);
    fHist_AntiDeuteron_TOF_MassSquare_pT->Fill(pT,massSq);
    fHist_AntiDeuteron_TOF_MassSquare_p->Fill(pTPC,massSq);

    // ITS
    fHist_AntiDeuteron_ITS_nCluster->Fill(Track->GetITSNcls());
    fHist_AntiDeuteron_ITS_nSigma_pT->Fill(pT,nSigmaITS);
    fHist_AntiDeuteron_ITS_nSigma_p->Fill(p,nSigmaITS);
    fHist_AntiDeuteron_ITS_dEdx_pT->Fill(pT,Track->GetITSsignal());
    fHist_AntiDeuteron_ITS_dEdx_p->Fill(p,Track->GetITSsignal());

    AntiDeuteronIsSelected = true;

    // save details of selected deuteron
    AntiDeuteronTrackArray->push_back(track);

  } // end of antideuteron loop





  // number of events containing selected ...
  if(AntiProtonIsSelected)			      fHist_Event_CutCounter->Fill(4);
  if(AntiDeuteronIsSelected)			      fHist_Event_CutCounter->Fill(5);
  if(AntiProtonIsSelected && AntiDeuteronIsSelected)  fHist_Event_CutCounter->Fill(6);


  int nAntiProtonsSelected	    = AntiProtonTrackArray->size();
  int nAntiDeuteronsSelected	    = AntiDeuteronTrackArray->size();
  int nAntiProtonAntiDeuteronPairs  = 0;





  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antiproton-antideuteron pairing loop ++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(AntiProtonIsSelected && AntiDeuteronIsSelected)
  {

    for(int Track1 = 0; Track1 < nAntiProtonsSelected; Track1++)	// particle pair loop (antiproton loop)
    { 

      AliAODTrack *AntiProtonTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(AntiProtonTrackArray->at(Track1)));

      double MomentumAntiProton[3];
      MomentumAntiProton[0] = AntiProtonTrack->Px();
      MomentumAntiProton[1] = AntiProtonTrack->Py();
      MomentumAntiProton[2] = AntiProtonTrack->Pz();

      TLorentzVector LorentzVectorAntiProton;
      LorentzVectorAntiProton.SetXYZM(MomentumAntiProton[0],MomentumAntiProton[1],MomentumAntiProton[2],ProtonMass);

      for(int Track2 = 0; Track2 < nAntiDeuteronsSelected; Track2++) // particle pair loop (antideuteron loop)
      {

	AliAODTrack *AntiDeuteronTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(AntiDeuteronTrackArray->at(Track2)));

	double MomentumAntiDeuteron[3];
	MomentumAntiDeuteron[0] = AntiDeuteronTrack->Px();
	MomentumAntiDeuteron[1] = AntiDeuteronTrack->Py();
	MomentumAntiDeuteron[2] = AntiDeuteronTrack->Pz();

	TLorentzVector LorentzVectorAntiDeuteron;
	LorentzVectorAntiDeuteron.SetXYZM(MomentumAntiDeuteron[0],MomentumAntiDeuteron[1],MomentumAntiDeuteron[2],DeuteronMass);

	fHist_AntiProtonAntiDeuteron_CutCounter->Fill(0);

	int AntiProtonID = AntiProtonTrack->GetID();
	int AntiDeuteronID = AntiDeuteronTrack->GetID();

	if((AntiProtonID == AntiDeuteronID) && DebugPairSelection)
	{
      
	  cout << "x-x-x-x-x-x-x-x-> same ID for proton and deuteron: " << AntiProtonID << endl;
	  cout << "AntiProton:\t\t\tAntiDeuteron:" << endl;
	  cout << MomentumAntiProton[0] << "\t\t\t" << MomentumAntiDeuteron[0] << endl;
	  cout << MomentumAntiProton[1] << "\t\t\t" << MomentumAntiDeuteron[1] << endl;
	  cout << MomentumAntiProton[2] << "\t\t\t" << MomentumAntiDeuteron[2] << endl;

	}

	if(AntiProtonID == AntiDeuteronID) continue;
	fHist_AntiProtonAntiDeuteron_CutCounter->Fill(1);

	// apply eta and phi cut
	bool RejectClosePairs = ClosePairRejection(bfield,*AntiProtonTrack,*AntiDeuteronTrack);
	if(RejectClosePairs) continue;
	fHist_AntiProtonAntiDeuteron_CutCounter->Fill(2);

	double angle = LorentzVectorAntiProton.Angle(LorentzVectorAntiDeuteron.Vect());

	// do pair cleaning (avoid auto-correlations)
	// apply pair cuts

	TLorentzVector LorentzVectorPair = LorentzVectorAntiProton + LorentzVectorAntiDeuteron;
	double RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorAntiProton,LorentzVectorAntiDeuteron);

	if(DebugPairSelection) cout << "angle: " << angle*TMath::RadToDeg() << "°\t kstar: " << RelativeMomentum << endl;
  
	nAntiProtonAntiDeuteronPairs++;
	fHist_AntiProtonAntiDeuteron_SED->Fill(RelativeMomentum);
	fHist_AntiProtonAntiDeuteron_AngleOfPairs->Fill(angle*TMath::RadToDeg(),RelativeMomentum);
	fHist_AntiProtonAntiDeuteron_Centrality->Fill(Centrality,RelativeMomentum);
	fHist_AntiProtonAntiDeuteron_VertexZ->Fill(PrimaryVertexZ,RelativeMomentum);

	if(RelativeMomentum <= kstar_max){

	  fHist_AntiProtonAntiDeuteron_pT->Fill(AntiProtonTrack->Pt(),AntiDeuteronTrack->Pt());
	  fHist_AntiProtonAntiDeuteron_Eta->Fill(AntiProtonTrack->Eta(),AntiDeuteronTrack->Eta());

	}


	double phi = AntiProtonTrack->Phi();

	if(DebugPairRotation){

	  cout << "---------------------------------" << endl;
	  cout << "original phi of proton: " << phi << " rad" << endl;
	  cout << "original k* of p-d pair: " << RelativeMomentum << endl;
	  cout << "original Px of proton: " << MomentumAntiProton[0] << endl;
	  cout << "original Py of proton: " << MomentumAntiProton[1] << endl;
	  cout << "original Pz of proton: " << MomentumAntiProton[2] << endl;
	  cout << "" << endl;

	} // end of DebugPairRotation

	// set up random generator
	TRandom3 *RandomGenerator = new TRandom3();
	RandomGenerator->SetSeed(0);

	// rotate pairs
	for(int iPair = 0; iPair < nRotatedPairs; iPair++)
	{

	  double RandomPhi = RandomGenerator->Uniform(0.0,2*TMath::Pi());
	  double phi_new = phi + RandomPhi;
	  AntiProtonTrack->SetPhi(phi_new);
      
	  double MomentumAntiProtonRotated[3];
	  MomentumAntiProtonRotated[0] = AntiProtonTrack->Px();
	  MomentumAntiProtonRotated[1] = AntiProtonTrack->Py();
	  MomentumAntiProtonRotated[2] = AntiProtonTrack->Pz();

	  TLorentzVector LorentzVectorAntiProtonRotated; 
	  LorentzVectorAntiProtonRotated.SetXYZM(MomentumAntiProtonRotated[0],MomentumAntiProtonRotated[1],MomentumAntiProtonRotated[2],ProtonMass);

	  TLorentzVector LorentzVectorPairRotated = LorentzVectorAntiProtonRotated + LorentzVectorAntiDeuteron;
	  double RelativeMomentumRotated = CalculateRelativeMomentum(LorentzVectorPairRotated,LorentzVectorAntiProtonRotated,LorentzVectorAntiDeuteron);

	  fHist_AntiProtonAntiDeuteron_RPD->Fill(RelativeMomentumRotated);

	  if(DebugPairRotation){

	    cout << "iPair: " << iPair << endl;
	    cout << "RandomPhi: " << RandomPhi << endl;
	    cout << "Phi_new: " << phi_new << endl;
	    cout << "Rotated k*: " << RelativeMomentumRotated << endl;
	    cout << "rotated Px of proton: " << MomentumAntiProtonRotated[0] << endl;
	    cout << "rotated Py of proton: " << MomentumAntiProtonRotated[1] << endl;
	    cout << "rotated Pz of proton: " << MomentumAntiProtonRotated[2] << endl;

	  } // end of DebugPairRotation

	} // end of pair rotation

      } // end of particle pair loop (deuteron loop)

    } // end of particle pair loop (proton loop)

  } // end of if statement: AntiProtonIsSelected && AntiDeuteronIsSelected


  if(DebugPairSelection)
  {

    if(!DebugEventSelection) cout << "" << endl;
    cout << "nAntiProtons selected:\t\t" << nAntiProtonsSelected << endl;
    cout << "nAntiDeuterons selected:\t\t" << nAntiDeuteronsSelected << endl;
    cout << "nPairs selected:\t\t" << nAntiProtonAntiDeuteronPairs << endl;

  } // end of DebugPairSelection

  fHist_AntiProtonAntiDeuteron_PairsPerEvent->Fill(nAntiProtonAntiDeuteronPairs);


  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ antiproton-antideuteron event mixing ++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(AntiProtonIsSelected && AntiDeuteronIsSelected) // fill antideuterons in event pool
  {

    int nMixedEvents  = EventPool_AntiDeuteron->GetCurrentNEvents();
    double ZvtxMin    = EventPool_AntiDeuteron->GetZvtxMin();
    double ZvtxMax    = EventPool_AntiDeuteron->GetZvtxMax();
    double MultMin    = EventPool_AntiDeuteron->GetMultMin();
    double MultMax    = EventPool_AntiDeuteron->GetMultMax();
    double PtMin      = EventPool_AntiDeuteron->GetPtMin();
    double PtMax      = EventPool_AntiDeuteron->GetPtMax();

    if(DebugEventMixing)
    {
    
      cout << "Number of events in buffer: " << nMixedEvents << endl;
      cout << "ZvtxMin: " << ZvtxMin << "\t\t ZvtxMax: " << ZvtxMax << endl;
      cout << "MultMin: " << MultMin << "\t\t MultMax: " << MultMax << endl;
      cout << "PtMin: " << PtMin << "\t\t PtMax: " << PtMax << endl;

    } // end of DebugEventMixing

    if(nMixedEvents == 0) fHist_AntiProtonAntiDeuteron_EventsForMixing->Fill(0);
    if(nMixedEvents > 0)  fHist_AntiProtonAntiDeuteron_EventsForMixing->Fill(1);

    for(int iMixedEvents = 0; iMixedEvents < nMixedEvents; iMixedEvents++){ // loop over mixed event

      TObjArray *MixedEventArray = (TObjArray*)EventPool_AntiDeuteron->GetEvent(iMixedEvents);
      if(!MixedEventArray) continue;

      AliAODTrackTiny *AntiDeuteronTrack = (AliAODTrackTiny*) MixedEventArray->At(0); 
      if(!AntiDeuteronTrack) continue;

      fHist_AntiProtonAntiDeuteron_UsedEventsInPool->Fill(AntiDeuteronTrack->GetCentrality(),AntiDeuteronTrack->GetPrimaryVertexZ());

      double MomentumAntiDeuteron[3];
      MomentumAntiDeuteron[0] = AntiDeuteronTrack->Px();
      MomentumAntiDeuteron[1] = AntiDeuteronTrack->Py();
      MomentumAntiDeuteron[2] = AntiDeuteronTrack->Pz();
   
      TLorentzVector LorentzVectorAntiDeuteron;
      LorentzVectorAntiDeuteron.SetXYZM(MomentumAntiDeuteron[0],MomentumAntiDeuteron[1],MomentumAntiDeuteron[2],DeuteronMass);

      for(int Track2 = 0; Track2 < nAntiProtonsSelected; Track2++) // particle pair loop (antiproton loop)
      {

	AliAODTrack *AntiProtonTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(AntiProtonTrackArray->at(Track2)));
	if(!AntiProtonTrack) continue;

	double MomentumAntiProton[3];
	MomentumAntiProton[0] = AntiProtonTrack->Px();
	MomentumAntiProton[1] = AntiProtonTrack->Py();
	MomentumAntiProton[2] = AntiProtonTrack->Pz();

	TLorentzVector LorentzVectorAntiProton;
	LorentzVectorAntiProton.SetXYZM(MomentumAntiProton[0],MomentumAntiProton[1],MomentumAntiProton[2],ProtonMass);

	TLorentzVector LorentzVectorPair = LorentzVectorAntiProton + LorentzVectorAntiDeuteron;
	double RelativeMomentum = CalculateRelativeMomentum(LorentzVectorPair,LorentzVectorAntiProton,LorentzVectorAntiDeuteron);

	if(TMath::IsNaN(RelativeMomentum)) continue;

	fHist_AntiProtonAntiDeuteron_MED->Fill(RelativeMomentum);

      } // loop over antiprotons

    } // loop over mixed events


    TObjArray *AntiDeuteronObjectArray = new TObjArray();
    AntiDeuteronObjectArray->SetOwner(true);

    int PickAntiDeuteronNumber = 0;

    if(nAntiDeuteronsSelected > 1) // randomize antideuteorn selection if more than one antideuteron is present in the current event
    {

      TRandom3 *RandomGenerator = new TRandom3();
      RandomGenerator->SetSeed(0);
      PickAntiDeuteronNumber = RandomGenerator->Integer(nAntiDeuteronsSelected); // produces a random integer between 0 and nAntiDeuteronsSelected-1

    }

    if(DebugRandomization) std::cout << "nAntiDeuteronsSelected: " << nAntiDeuteronsSelected << "\t PickAntiDeuteronNumber: " << PickAntiDeuteronNumber << std::endl;

    AliAODTrack *AntiDeuteronTrack = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(AntiDeuteronTrackArray->at(PickAntiDeuteronNumber))); // array starts at 0
    if(!AntiDeuteronTrack) std::cout << "AliAnalysisTask_pdLd: No AntiDeuteronTrack received from the Eventpool" << std::endl;

    AliAODTrackTiny *TinyAntiDeuteronTrack = new AliAODTrackTiny();
    TinyAntiDeuteronTrack->InitFromTrack(AntiDeuteronTrack,Centrality,PrimaryVertexZ);

    AntiDeuteronObjectArray->Add(TinyAntiDeuteronTrack);

    EventPool_AntiDeuteron->UpdatePool(AntiDeuteronObjectArray);

  } // end of antiproton and antideuteron is selected






  


  fHist_ProtonDeuteron_PairMultiplicity->Fill(DeuteronTrackArray->size(),ProtonTrackArray->size());
  fHist_AntiProtonAntiDeuteron_PairMultiplicity->Fill(AntiDeuteronTrackArray->size(),AntiProtonTrackArray->size());


  PostData(1,fHistList_Event);
  PostData(2,fHistList_Proton);
  PostData(3,fHistList_Deuteron);
  PostData(4,fHistList_ProtonDeuteron);
  PostData(5,fHistList_AntiProton);
  PostData(6,fHistList_AntiDeuteron);
  PostData(7,fHistList_AntiProtonAntiDeuteron);

} // end of UserExec















void AliAnalysisTask_pdLd::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
double AliAnalysisTask_pdLd::CalculateBetaTOF(AliAODTrack &track)
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
double AliAnalysisTask_pdLd::CalculateMassSquareTOF(AliAODTrack &track)
{

  double p = track.P();
  double beta = CalculateBetaTOF(track);
  double mass2 = -999.0;

  if(beta > 0.0){

    mass2 = (1/(beta*beta)-1) * (p*p);

  }

  return mass2;

} // end of CalculateMassSquareTOF





// calculate the relative momentum k*
double AliAnalysisTask_pdLd::CalculateRelativeMomentum(TLorentzVector &Pair, TLorentzVector &Part1, TLorentzVector &Part2)
{

  double beta = Pair.Beta();
  double betaX = beta * cos(Pair.Phi()) * sin(Pair.Theta());
  double betaY = beta * sin(Pair.Phi()) * sin(Pair.Theta());
  double betaZ = beta * cos(Pair.Theta());

  TLorentzVector Part1CMS = Part1;
  TLorentzVector Part2CMS = Part2;

  Part1CMS.Boost(-betaX,-betaY,-betaZ);
  Part2CMS.Boost(-betaX,-betaY,-betaZ);

  TLorentzVector RelativeMomentum = Part1CMS - Part2CMS;

  return 0.5 * RelativeMomentum.P();

} // end of CalculateRelativeMomentum


bool AliAnalysisTask_pdLd::ClosePairRejection(double bfield, AliAODTrack &track1, AliAODTrack &track2)
{

  bool reject = false;

  double minimum = 0.017;
  double minimum_sq = minimum * minimum;
  double radii[9] = {85,105,125,145,165,185,205,225,245}; // cm

  double eta1 = track1.Eta();
  double eta2 = track2.Eta();

  for(int radius = 0; radius < 9; radius++)
  {

    double phi1 = RecalculatePhi(bfield,track1,radii[radius]);
    double phi2 = RecalculatePhi(bfield,track2,radii[radius]);

    double delta_eta = eta1 - eta2;
    double delta_phi = phi1 - phi2;

    double delta_eta_sq = delta_eta*delta_eta;
    double delta_phi_sq = delta_phi*delta_phi;

    if((delta_eta_sq + delta_phi_sq) < (minimum_sq + minimum_sq)) reject = true;

  }

  return reject;

} // end of ClosePairRejection


double AliAnalysisTask_pdLd::RecalculatePhi(double bfield, AliAODTrack &track, double radius)
{

  double phi0	= track.Phi();
  double charge = track.Charge();
  double pT	= track.Pt();

  double argument = ((0.1 * bfield) * charge * (radius * 0.01) * 0.3)/(2 * pT);
  double phi = phi0 - TMath::ASin(TMath::Abs(argument));

  return phi;

} // end of RecalculatePhi




double AliAnalysisTask_pdLd::CalculateDeuteronSigmaMassSquareTOF(double pT, double massSq, bool isMatter)
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
bool AliAnalysisTask_pdLd::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter)
{

  bool PassedParticleCuts = false;

  // define proton and antiproton track cuts
  double Proton_pT_min = 0.5;
  double Proton_pT_max = 4.0;
  double Proton_eta_min = -0.8;
  double Proton_eta_max = +0.8;
  double Proton_DCAxy_max = 0.1; // cm
  double Proton_DCAz_max = 0.2; // cm

  double Proton_TPC_RatioRowsCluster_min = 0.83;
  double Proton_TPC_nSigma_max = 3.0;
  double Proton_TOF_nSigma_max = 3.0;
  double Proton_TOF_nSigma_max_low_p = 7.0;
  int Proton_TPC_nCluster_min = 80;
  int Proton_TPC_nCrossedRows_min = 70;
  int Proton_TPC_nSharedCluster_max = 0;
  int Proton_FilterBit = 7; // FilterBit 7 = 128
  int Proton_ITS_nCluster_min = 2;
  double Proton_ITS_nSigma_max = 60;

  bool UseITS = true;
  double Proton_TPC_Threshold = 0.70;


  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(1);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(1);

  // apple TPC sigma cut
  double nSigmaTPC = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kProton);
  if(TMath::Abs(nSigmaTPC) > Proton_TPC_nSigma_max) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(2);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(2);

  // get DCA infotmation
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  double DCAxy = xv[0];
  double DCAz = xv[1];

  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Proton_DCAxy_max) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(3);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(3);

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Proton_DCAz_max) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(4);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(4);

  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  double p    = Track.P();
  double pT   = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  // check if TOF information is available above threshold
  if((pTPC >= Proton_TPC_Threshold) && (!TOFisOK)) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(5);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(5);

  double nSigmaTOF = -999.0;
  if(TOFisOK) nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kProton);

  // apply TOF low p cut
  if(TOFisOK)
  {

    if((pTPC < Proton_TPC_Threshold) && (TMath::Abs(nSigmaTOF) > Proton_TOF_nSigma_max_low_p)) return PassedParticleCuts;

  }
  if(isMatter)	fHist_Proton_CutCounter->Fill(6);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(6);

  // apply TOF high p cut
  if(TOFisOK)
  {

    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kProton);

    if((pTPC >= Proton_TPC_Threshold) && (TMath::Abs(nSigmaTOF) > Proton_TOF_nSigma_max)) return PassedParticleCuts;

  }
  if(isMatter)	fHist_Proton_CutCounter->Fill(7);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(7);

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
  if(isMatter)	fHist_Proton_CutCounter->Fill(8);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(8);
*/

  // apply FilterBit cut
  if(!Track.TestFilterBit(Proton_FilterBit)) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(9);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(9);

  // apply pT cut
  if(pT < Proton_pT_min || pT > Proton_pT_max) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(10);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(10);

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(11);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(11);

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(eta < Proton_eta_min || eta > Proton_eta_max) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(12);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(12);

  // apply cluster cut for TPC
  int nClusterTPC = Track.GetNcls(1);
  if(nClusterTPC < Proton_TPC_nCluster_min) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(13);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(13);

  // apply crossed rows cut for TPC
  int nCrossedRowsTPC = Track.GetTPCCrossedRows();
  if(nCrossedRowsTPC < Proton_TPC_nCrossedRows_min) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(14);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(14);

  // apply zero shared cluster cut for TPC
  int nSharedClusterTPC = Track.GetTPCnclsS();
  if(nSharedClusterTPC > Proton_TPC_nSharedCluster_max) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(15);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(15);

  // apply findable cluster cut for TPC
  int nClusterTPCfindable = Track.GetTPCNclsF();
  double RatioRowsFindableClusterTPC = -999.0;
  if(nClusterTPCfindable > 0) RatioRowsFindableClusterTPC = ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);
  if(RatioRowsFindableClusterTPC < Proton_TPC_RatioRowsCluster_min) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(16);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(16);


  // apply ITS cluster cut
  double nClusterITS = Track.GetITSNcls();
  if((UseITS) && (nClusterITS <= Proton_ITS_nCluster_min)) return PassedParticleCuts;
  if(isMatter)	fHist_Proton_CutCounter->Fill(17);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(17);

  // apply ITS dEdx cut below proton band
  if((UseITS))
  {

    TF1 *fdEdxProtonITS = new TF1("fdEdxProtonITS","0.7*[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.);
    fdEdxProtonITS->SetParameters(2.36861e-07,-55831.1,-238672,9.55834,17081,1,0.93827208816);

    if(Track.GetITSsignal() < fdEdxProtonITS->Eval(p)) return PassedParticleCuts;

  }

  // apply ITS dEdx cut above proton band
  if((UseITS))
  {

    TF1 *fdEdxProtonITS = new TF1("fdEdxProtonITS","1.3*[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.);
    fdEdxProtonITS->SetParameters(2.36861e-07,-55831.1,-238672,9.55834,17081,1,0.93827208816);

    if(Track.GetITSsignal() > fdEdxProtonITS->Eval(p)) return PassedParticleCuts;

  }
  if(isMatter)	fHist_Proton_CutCounter->Fill(18);
  if(!isMatter) fHist_AntiProton_CutCounter->Fill(18);


  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckProtonCuts


// apply track cuts for deuterons and antideuterons
bool AliAnalysisTask_pdLd::CheckDeuteronCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter)
{

  bool PassedParticleCuts = false;

  // define deuteron and antideuteron track cuts
  double Deuteron_pT_min = 0.5;
  double Deuteron_pT_max = 1.4;
  double Deuteron_eta_min = -0.8;
  double Deuteron_eta_max = +0.8;
  double Deuteron_DCAxy_max = 0.1; // cm
  double Deuteron_DCAz_max = 0.2; // cm

  double Deuteron_TPC_RatioRowsCluster_min = 0.83;
  double Deuteron_TPC_nSigma_max = 3.0;
  double Deuteron_TOF_nSigma_max = 3.0;
  double Deuteron_TOF_nSigma_max_low_p = 7.0;
  int Deuteron_TPC_nCluster_min = 80;
  int Deuteron_TPC_nCrossedRows_min = 70;
  int Deuteron_TPC_nSharedCluster_max = 0;
  int Deuteron_FilterBit = 7; // FilterBit 7 = 128
  int Deuteron_ITS_nCluster_min = 2;
  double Deuteron_ITS_nSigma_max = 60;


  double Deuteron_TPC_Threshold = 1.0;

  double Deuteron_TOF_MassSquare_nSigma_max = 3.0;

  bool UseBetaTOF = false;
  bool UseMassSquareTOF = true;
  bool UseITS = true;


  // check if TPC information is available
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTPC,&Track);
  if(!(statusTPC == AliPIDResponse::kDetPidOk)) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(1);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(1);

  // apply TPC nSigma cut
  double nSigmaTPC = fPIDResponse.NumberOfSigmasTPC(&Track,AliPID::kDeuteron);
  if(TMath::Abs(nSigmaTPC) > Deuteron_TPC_nSigma_max) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(2);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(2);

  // get DCA information
  float xv[2];
  float yv[3];
  Track.GetImpactParameters(xv,yv);
  double DCAxy = xv[0];
  double DCAz = xv[1];
  
  // apply DCAxy cut
  if(TMath::Abs(DCAxy) > Deuteron_DCAxy_max) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(3);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(3);

  // apply DCAz cut
  if(TMath::Abs(DCAz) > Deuteron_DCAz_max) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(4);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(4);

  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse.CheckPIDStatus(AliPIDResponse::kTOF,&Track);
  bool TOFisOK = false;
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

  double p = Track.P();
  double pT = Track.Pt();
  double pTPC = Track.GetTPCmomentum();

  // check TOF status above threshold
  if((pTPC >= Deuteron_TPC_Threshold) && (!TOFisOK)) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(5);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(5);

  double nSigmaTOF = -999.0;
    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);

  // apply TOF low p cut
  if(TOFisOK)
  {

    if((pTPC < Deuteron_TPC_Threshold) && (TMath::Abs(nSigmaTOF) > Deuteron_TOF_nSigma_max_low_p)) return PassedParticleCuts;

  }
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(6);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(6);


  if((UseBetaTOF) && (pTPC >= Deuteron_TPC_Threshold))
  {

    // apply TPC and TOF sigma cut

    nSigmaTOF = fPIDResponse.NumberOfSigmasTOF(&Track,AliPID::kDeuteron);
    if((TMath::Abs(nSigmaTPC) > Deuteron_TPC_nSigma_max) || (TMath::Abs(nSigmaTOF) > Deuteron_TOF_nSigma_max)) return PassedParticleCuts;

    if(isMatter)  fHist_Deuteron_CutCounter->Fill(6);
    if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(6);


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

  if(isMatter)	fHist_Deuteron_CutCounter->Fill(7);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(7);

  } // end of UseBetaTOF


  if(UseMassSquareTOF)
  {
  
    if(pTPC >= Deuteron_TPC_Threshold)
    {

    double massSq = CalculateMassSquareTOF(Track);
    double nSigmaTOFmSq = CalculateDeuteronSigmaMassSquareTOF(pT,massSq,isMatter);

    if(TMath::Abs(nSigmaTOFmSq) > Deuteron_TOF_MassSquare_nSigma_max) return PassedParticleCuts;

    } // end of UseMassSquareTOF

    if(isMatter)  fHist_Deuteron_CutCounter->Fill(9);
    if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(9);

  }

  // apply FilterBit cut
  if(!Track.TestFilterBit(Deuteron_FilterBit)) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(10);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(10);

  // apply pT cut
  if(pT < Deuteron_pT_min || pT > Deuteron_pT_max) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(11);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(11);

  // apply charge cut
  int charge = Track.Charge();
  if(charge < 1 && isMatter)   return PassedParticleCuts;
  if(charge > -1 && !isMatter) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(12);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(12);

  // apply pseudo-rapidity cut
  double eta = Track.Eta();
  if(eta < Deuteron_eta_min || eta > Deuteron_eta_max) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(13);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(13);

  // apply cluster cut for TPC
  int nClusterTPC = Track.GetNcls(1);
  if(nClusterTPC < Deuteron_TPC_nCluster_min) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(14);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(14);

  // apply crossed rows cut for TPC
  int nCrossedRowsTPC = Track.GetTPCCrossedRows();
  if(nCrossedRowsTPC < Deuteron_TPC_nCrossedRows_min) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(15);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(15);

  // apply zero shared cluster cut for TPC
  int nSharedClusterTPC = Track.GetTPCnclsS();
  if(nSharedClusterTPC > Deuteron_TPC_nSharedCluster_max) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(16);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(16);

  // apply findable cluster cut for TPC
  int nClusterTPCfindable = Track.GetTPCNclsF();
  double RatioRowsFindableClusterTPC = -999.0;
  if(nClusterTPCfindable > 0) RatioRowsFindableClusterTPC = ((double)nCrossedRowsTPC / (double)nClusterTPCfindable);
  if(RatioRowsFindableClusterTPC < Deuteron_TPC_RatioRowsCluster_min) return PassedParticleCuts;
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(17);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(17);


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
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(21);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(21);


  // apply ITS dEdx cut
  if(UseITS)
  {

    TF1 *fdEdxDeuteronITS = new TF1("fdEdxDeuteronITS","0.7*[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.1,6.);
    fdEdxDeuteronITS->SetParameters(7.41722e-06,-55831.1,-238672,11249.3,19828.9,1,1.8756129425);

    if(Track.GetITSsignal() < fdEdxDeuteronITS->Eval(p)) return PassedParticleCuts;

  }
  if(isMatter)	fHist_Deuteron_CutCounter->Fill(20);
  if(!isMatter) fHist_AntiDeuteron_CutCounter->Fill(20);



  PassedParticleCuts = true;
  return PassedParticleCuts;

} // end of CheckDeuteronCuts



