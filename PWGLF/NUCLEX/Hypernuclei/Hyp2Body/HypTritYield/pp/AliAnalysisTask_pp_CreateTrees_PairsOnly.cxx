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

#include "AliAnalysisTask_pp_CreateTrees_PairsOnly.h"

using namespace std;
ClassImp(AliAnalysisTask_pp_CreateTrees_PairsOnly) 






AliAnalysisTask_pp_CreateTrees_PairsOnly::AliAnalysisTask_pp_CreateTrees_PairsOnly() : AliAnalysisTaskSE(),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(0),
  fUseOpenCuts(0),
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
  fProton_ID(0),
  fProton_Event_Multiplicity(0),
  fProton_Event_Identifier(0),
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
  fAntiProton_ID(0),
  fAntiProton_Event_Multiplicity(0),
  fAntiProton_Event_Identifier(0),
  fHistoList(0),
  h_Proton_TOF_m2_NoTOFcut(0),
  h_AntiProton_TOF_m2_NoTOFcut(0),
  h_Proton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut(0)
{


}



AliAnalysisTask_pp_CreateTrees_PairsOnly::AliAnalysisTask_pp_CreateTrees_PairsOnly(const char *name,int CollisionSystem, bool UseOpenCuts) : AliAnalysisTaskSE(name),
  fAODEvent(0),
  fAODHandler(0),
  fHeader(0),
  fPIDResponse(0),
  fCollisionSystem(CollisionSystem),
  fUseOpenCuts(UseOpenCuts),
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
  fProton_ID(0),
  fProton_Event_Multiplicity(0),
  fProton_Event_Identifier(0),
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
  fAntiProton_ID(0),
  fAntiProton_Event_Multiplicity(0),
  fAntiProton_Event_Identifier(0),
  fHistoList(0),
  h_Proton_TOF_m2_NoTOFcut(0),
  h_AntiProton_TOF_m2_NoTOFcut(0),
  h_Proton_ITS_dEdx_NoTOFcutNoITScut(0),
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut(0)
{

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TList::Class());

}

  
AliAnalysisTask_pp_CreateTrees_PairsOnly::~AliAnalysisTask_pp_CreateTrees_PairsOnly()
{

  if(fSaveTree_Proton)
    {
      delete fSaveTree_Proton;
    }

  if(fSaveTree_AntiProton)
    {
      delete fSaveTree_AntiProton;
    }

  if(fHistoList)
    {
      delete fHistoList;
    }

}







void AliAnalysisTask_pp_CreateTrees_PairsOnly::UserCreateOutputObjects()
{

  fHistoList = new TList();
  fHistoList->SetOwner();

  h_Proton_TOF_m2_NoTOFcut = new TH2F("h_Proton_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (protons)",240,0.0,6.0,500,0.0,10.0);
  h_Proton_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_Proton_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_Proton_TOF_m2_NoTOFcut);

  h_AntiProton_TOF_m2_NoTOFcut = new TH2F("h_AntiProton_TOF_m2_NoTOFcut","TOF #it{m}^{2} without TOF cut (antiprotons)",240,0.0,6.0,500,0.0,10.0);
  h_AntiProton_TOF_m2_NoTOFcut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h_AntiProton_TOF_m2_NoTOFcut->GetYaxis()->SetTitle("#it{m}^{2} (GeV/#it{c}^{2})^{2}");
  fHistoList->Add(h_AntiProton_TOF_m2_NoTOFcut);


  h_Proton_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_Proton_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (protons)",400,0.0,10.0,500,0.0,500.0);
  h_Proton_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_Proton_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_Proton_ITS_dEdx_NoTOFcutNoITScut);

  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut = new TH2F("h_AntiProton_ITS_dEdx_NoTOFcutNoITScut","ITS d#it{E}/d#it{x} without TOF and ITS cut (antiprotons)",400,0.0,10.0,500,0.0,500.0);
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  h_AntiProton_ITS_dEdx_NoTOFcutNoITScut->GetYaxis()->SetTitle("d#it{E}/d#it{x} (abs.)");
  fHistoList->Add(h_AntiProton_ITS_dEdx_NoTOFcutNoITScut);


  fSaveTree_Proton = new TTree("fSaveTree_Proton","fSaveTree_Proton");
  fSaveTree_Proton->Branch("Proton_px",&fProton_px,"Proton_px/F");
  fSaveTree_Proton->Branch("Proton_py",&fProton_py,"Proton_py/F");
  fSaveTree_Proton->Branch("Proton_pz",&fProton_pz,"Proton_pz/F");
  fSaveTree_Proton->Branch("Proton_pTPC",&fProton_pTPC,"Proton_pTPC/F");
  fSaveTree_Proton->Branch("Proton_Eta",&fProton_Eta,"Proton_Eta/F");
  fSaveTree_Proton->Branch("Proton_Phi",&fProton_Phi,"Proton_Phi/F");
  fSaveTree_Proton->Branch("Proton_TPC_Chi2",&fProton_TPC_Chi2,"Proton_TPC_Chi2/F");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx",&fProton_TPC_dEdx,"Proton_TPC_dEdx/F");
  fSaveTree_Proton->Branch("Proton_TPC_dEdx_nSigma",&fProton_TPC_dEdx_nSigma,"Proton_TPC_dEdx_nSigma/F");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2",&fProton_TOF_Mass2,"Proton_TOF_Mass2/F");
  fSaveTree_Proton->Branch("Proton_TOF_Mass2_nSigma",&fProton_TOF_Mass2_nSigma,"Proton_TOF_Mass2_nSigma/F");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx",&fProton_ITS_dEdx,"Proton_ITS_dEdx/F");
  fSaveTree_Proton->Branch("Proton_ITS_dEdx_nSigma",&fProton_ITS_dEdx_nSigma,"Proton_ITS_dEdx_nSigma/F");
  fSaveTree_Proton->Branch("Proton_DCAxy",&fProton_DCAxy,"Proton_DCAxy/F");
  fSaveTree_Proton->Branch("Proton_DCAz",&fProton_DCAz,"Proton_DCAz/F");
  fSaveTree_Proton->Branch("Proton_Event_Centrality",&fProton_Event_Centrality,"Proton_Event_Centrality/F");
  fSaveTree_Proton->Branch("Proton_Event_PrimaryVertexZ",&fProton_Event_PrimaryVertexZ,"Proton_Event_PrimaryVertexZ/F");
  fSaveTree_Proton->Branch("Proton_Event_BField",&fProton_Event_BField,"Proton_Event_BField/F");
  fSaveTree_Proton->Branch("Proton_TPC_nCrossedRows",&fProton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fSaveTree_Proton->Branch("Proton_TPC_nSharedCluster",&fProton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fSaveTree_Proton->Branch("Proton_TPC_nFindableCluster",&fProton_TPC_nFindableCluster,"Proton_TPC_nFindableCluster/s");
  fSaveTree_Proton->Branch("Proton_TPC_nCluster",&fProton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fSaveTree_Proton->Branch("Proton_ITS_nCluster",&fProton_ITS_nCluster,"Proton_ITS_nCluster/s");
  fSaveTree_Proton->Branch("Proton_ID",&fProton_ID,"Proton_ID/i");
  fSaveTree_Proton->Branch("Proton_Event_Multiplicity",&fProton_Event_Multiplicity,"Proton_Event_Multiplicity/i");
  fSaveTree_Proton->Branch("Proton_Event_Identifier",&fProton_Event_Identifier,"Proton_Event_Identifier/l");




  fSaveTree_AntiProton = new TTree("fSaveTree_AntiProton","fSaveTree_AntiProton");
  fSaveTree_AntiProton->Branch("AntiProton_px",&fAntiProton_px,"AntiProton_px/F");
  fSaveTree_AntiProton->Branch("AntiProton_py",&fAntiProton_py,"AntiProton_py/F");
  fSaveTree_AntiProton->Branch("AntiProton_pz",&fAntiProton_pz,"AntiProton_pz/F");
  fSaveTree_AntiProton->Branch("AntiProton_pTPC",&fAntiProton_pTPC,"AntiProton_pTPC/F");
  fSaveTree_AntiProton->Branch("AntiProton_Eta",&fAntiProton_Eta,"AntiProton_Eta/F");
  fSaveTree_AntiProton->Branch("AntiProton_Phi",&fAntiProton_Phi,"AntiProton_Phi/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_Chi2",&fAntiProton_TPC_Chi2,"AntiProton_TPC_Chi2/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx",&fAntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_dEdx_nSigma",&fAntiProton_TPC_dEdx_nSigma,"AntiProton_TPC_dEdx_nSigma/F");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2",&fAntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/F");
  fSaveTree_AntiProton->Branch("AntiProton_TOF_Mass2_nSigma",&fAntiProton_TOF_Mass2_nSigma,"AntiProton_TOF_Mass2_nSigma/F");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx",&fAntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/F");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_dEdx_nSigma",&fAntiProton_ITS_dEdx_nSigma,"AntiProton_ITS_dEdx_nSigma/F");
  fSaveTree_AntiProton->Branch("AntiProton_DCAxy",&fAntiProton_DCAxy,"AntiProton_DCAxy/F");
  fSaveTree_AntiProton->Branch("AntiProton_DCAz",&fAntiProton_DCAz,"AntiProton_DCAz/F");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Centrality",&fAntiProton_Event_Centrality,"AntiProton_Event_Centrality/F");
  fSaveTree_AntiProton->Branch("AntiProton_Event_PrimaryVertexZ",&fAntiProton_Event_PrimaryVertexZ,"AntiProton_Event_PrimaryVertexZ/F");
  fSaveTree_AntiProton->Branch("AntiProton_Event_BField",&fAntiProton_Event_BField,"AntiProton_Event_BField/F");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&fAntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&fAntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nFindableCluster",&fAntiProton_TPC_nFindableCluster,"AntiProton_TPC_nFindableCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_TPC_nCluster",&fAntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_ITS_nCluster",&fAntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/s");
  fSaveTree_AntiProton->Branch("AntiProton_ID",&fAntiProton_ID,"AntiProton_ID/i");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Multiplicity",&fAntiProton_Event_Multiplicity,"AntiProton_Event_Multiplicity/i");
  fSaveTree_AntiProton->Branch("AntiProton_Event_Identifier",&fAntiProton_Event_Identifier,"AntiProton_Event_Identifier/l");





  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_AntiProton);
  PostData(3,fHistoList);



} // end of UserCreateOutputObjects









void AliAnalysisTask_pp_CreateTrees_PairsOnly::UserExec(Option_t*)
{

//  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAODEvent)::Fatal("AliAnalysisTask_pp_CreateTrees_PairsOnly::UserExec","No AOD event found!");

  fHeader = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  if(!fHeader)::Fatal("AliAnalysisTask_pp_CreateTrees_PairsOnly::UserExec","No Header found!");

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse)::Fatal("AliAnalysisTask_pp_CreateTrees_PairsOnly::UserExec","No PIDResponse found!");


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
  if(!PrimaryVertex)::Warning("AliAnalsisTask_pp_CreateTrees_PairsOnlyd::UserExec","No AliAODVertex object found!");
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

  // EventID -> https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventInfo
  unsigned long EventID	= (unsigned long)BunchCrossNumber + ((unsigned long)OrbitNumber*3564) + ((unsigned long)PeriodNumber*16777215*3564);





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
    cout << "z-position of primary vertex:\t" << PrimaryVertexZ << " cm" << endl;

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
  unsigned short    Proton_ITS_nCluster;
  unsigned int      Proton_ID;
  unsigned long     Proton_Event_Identifier;


  TTree *fTempTree_Proton = new TTree("fTempTree_Proton","fTempTree_Proton");
  fTempTree_Proton->Branch("Proton_px",&Proton_px,"Proton_px/F");
  fTempTree_Proton->Branch("Proton_py",&Proton_py,"Proton_py/F");
  fTempTree_Proton->Branch("Proton_pz",&Proton_pz,"Proton_pz/F");
  fTempTree_Proton->Branch("Proton_pTPC",&Proton_pTPC,"Proton_pTPC/F");
  fTempTree_Proton->Branch("Proton_Eta",&Proton_Eta,"Proton_Eta/F");
  fTempTree_Proton->Branch("Proton_Phi",&Proton_Phi,"Proton_Phi/F");
  fTempTree_Proton->Branch("Proton_TPC_Chi2",&Proton_TPC_Chi2,"Proton_TPC_Chi2/F");
  fTempTree_Proton->Branch("Proton_TPC_dEdx",&Proton_TPC_dEdx,"Proton_TPC_dEdx/F");
  fTempTree_Proton->Branch("Proton_TPC_dEdx_nSigma",&Proton_TPC_dEdx_nSigma,"Proton_TPC_dEdx_nSigma/F");
  fTempTree_Proton->Branch("Proton_TOF_Mass2",&Proton_TOF_Mass2,"Proton_TOF_Mass2/F");
  fTempTree_Proton->Branch("Proton_TOF_Mass2_nSigma",&Proton_TOF_Mass2_nSigma,"Proton_TOF_Mass2_nSigma/F");
  fTempTree_Proton->Branch("Proton_ITS_dEdx",&Proton_ITS_dEdx,"Proton_ITS_dEdx/F");
  fTempTree_Proton->Branch("Proton_ITS_dEdx_nSigma",&Proton_ITS_dEdx_nSigma,"Proton_ITS_dEdx_nSigma/F");
  fTempTree_Proton->Branch("Proton_DCAxy",&Proton_DCAxy,"Proton_DCAxy/F");
  fTempTree_Proton->Branch("Proton_DCAz",&Proton_DCAz,"Proton_DCAz/F");
  fTempTree_Proton->Branch("Proton_TPC_nCrossedRows",&Proton_TPC_nCrossedRows,"Proton_TPC_nCrossedRows/s");
  fTempTree_Proton->Branch("Proton_TPC_nSharedCluster",&Proton_TPC_nSharedCluster,"Proton_TPC_nSharedCluster/s");
  fTempTree_Proton->Branch("Proton_TPC_nFindableCluster",&Proton_TPC_nFindableCluster,"Proton_TPC_nFindableCluster/s");
  fTempTree_Proton->Branch("Proton_TPC_nCluster",&Proton_TPC_nCluster,"Proton_TPC_nCluster/s");
  fTempTree_Proton->Branch("Proton_ITS_nCluster",&Proton_ITS_nCluster,"Proton_ITS_nCluster/s");
  fTempTree_Proton->Branch("Proton_ID",&Proton_ID,"Proton_ID/i");
  fTempTree_Proton->Branch("Proton_Event_Identifier",&Proton_Event_Identifier,"Proton_Event_Identifier/l");

  unsigned short nProtonsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// proton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pp_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;

    // apply proton cuts
    bool PassedProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,true,RunNumber,fUseOpenCuts);
    if(!PassedProtonCuts) continue;
  
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
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,1,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,1,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();
      if(ITS_nCluster < 0) ITS_nCluster = 0;

    }


    Proton_px			    = Track->Px();
    Proton_py			    = Track->Py();
    Proton_pz			    = Track->Pz();
    Proton_pTPC			    = Track->GetTPCmomentum();
    Proton_Eta			    = Track->Eta();
    Proton_Phi			    = Track->Phi();
    Proton_TPC_Chi2		    = Track->GetTPCchi2();
    Proton_TPC_dEdx		    = Track->GetTPCsignal();
    Proton_TPC_dEdx_nSigma	    = (float)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    Proton_TOF_Mass2		    = (float)TOF_m2;
    Proton_TOF_Mass2_nSigma	    = (float)TOF_m2_nSigma;
    Proton_ITS_dEdx		    = (float)ITS_dEdx;
    Proton_ITS_dEdx_nSigma	    = (float)ITS_dEdx_nSigma;
    Proton_DCAxy		    = DCAxy;
    Proton_DCAz			    = DCAz;
    Proton_TPC_nCrossedRows	    = Track->GetTPCCrossedRows();
    Proton_TPC_nSharedCluster	    = Track->GetTPCnclsS();
    Proton_TPC_nFindableCluster	    = Track->GetTPCNclsF();
    Proton_TPC_nCluster		    = Track->GetTPCNcls();
    Proton_ITS_nCluster		    = (unsigned short)ITS_nCluster;
    Proton_ID			    = track;
    Proton_Event_Identifier	    = EventID;
 
    fTempTree_Proton->Fill();
    nProtonsSelected++;

  } // end of proton loop




  if(nProtonsSelected > 1){

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

      fProton_Event_Multiplicity    = Multiplicity;
      fProton_Event_Centrality	    = Centrality;
      fProton_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fProton_Event_BField	    = BField;
      

      fSaveTree_Proton->Fill();

    } // end of loop (copy protons)


  }

  fTempTree_Proton->Delete();








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
  unsigned short    AntiProton_ITS_nCluster;
  unsigned int      AntiProton_ID;
  unsigned long     AntiProton_Event_Identifier;

  fTempTree_AntiProton = new TTree("fTempTree_AntiProton","fTempTree_AntiProton");
  fTempTree_AntiProton->Branch("AntiProton_px",&AntiProton_px,"AntiProton_px/F");
  fTempTree_AntiProton->Branch("AntiProton_py",&AntiProton_py,"AntiProton_py/F");
  fTempTree_AntiProton->Branch("AntiProton_pz",&AntiProton_pz,"AntiProton_pz/F");
  fTempTree_AntiProton->Branch("AntiProton_pTPC",&AntiProton_pTPC,"AntiProton_pTPC/F");
  fTempTree_AntiProton->Branch("AntiProton_Eta",&AntiProton_Eta,"AntiProton_Eta/F");
  fTempTree_AntiProton->Branch("AntiProton_Phi",&AntiProton_Phi,"AntiProton_Phi/F");
  fTempTree_AntiProton->Branch("AntiProton_TPC_Chi2",&AntiProton_TPC_Chi2,"AntiProton_TPC_Chi2/F");
  fTempTree_AntiProton->Branch("AntiProton_TPC_dEdx",&AntiProton_TPC_dEdx,"AntiProton_TPC_dEdx/F");
  fTempTree_AntiProton->Branch("AntiProton_TPC_dEdx_nSigma",&AntiProton_TPC_dEdx_nSigma,"AntiProton_TPC_dEdx_nSigma/F");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Mass2",&AntiProton_TOF_Mass2,"AntiProton_TOF_Mass2/F");
  fTempTree_AntiProton->Branch("AntiProton_TOF_Mass2_nSigma",&AntiProton_TOF_Mass2_nSigma,"AntiProton_TOF_Mass2_nSigma/F");
  fTempTree_AntiProton->Branch("AntiProton_ITS_dEdx",&AntiProton_ITS_dEdx,"AntiProton_ITS_dEdx/F");
  fTempTree_AntiProton->Branch("AntiProton_ITS_dEdx_nSigma",&AntiProton_ITS_dEdx_nSigma,"AntiProton_ITS_dEdx_nSigma/F");
  fTempTree_AntiProton->Branch("AntiProton_DCAxy",&AntiProton_DCAxy,"AntiProton_DCAxy/F");
  fTempTree_AntiProton->Branch("AntiProton_DCAz",&AntiProton_DCAz,"AntiProton_DCAz/F");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nCrossedRows",&AntiProton_TPC_nCrossedRows,"AntiProton_TPC_nCrossedRows/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nSharedCluster",&AntiProton_TPC_nSharedCluster,"AntiProton_TPC_nSharedCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nFindableCluster",&AntiProton_TPC_nFindableCluster,"AntiProton_TPC_nFindableCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_TPC_nCluster",&AntiProton_TPC_nCluster,"AntiProton_TPC_nCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_ITS_nCluster",&AntiProton_ITS_nCluster,"AntiProton_ITS_nCluster/s");
  fTempTree_AntiProton->Branch("AntiProton_ID",&AntiProton_ID,"AntiProton_ID/i");
  fTempTree_AntiProton->Branch("AntiProton_Event_Identifier",&AntiProton_Event_Identifier,"AntiProton_Event_Identifier/l");

  unsigned short nAntiProtonsSelected = 0;

  for(int track = 0; track < nTracks; track++)	// antiproton loop
  { 

    AliAODTrack *Track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(track));
    if(!Track)::Warning("AliAnalysisTask_pp_CreateTrees_PairsOnly::UserExec","No AliAODTrack found");
    if(!Track) continue;
  

    // apply antiproton cuts
    bool PassedAntiProtonCuts = CheckProtonCuts(*Track,*fPIDResponse,false,RunNumber,fUseOpenCuts);
    if(!PassedAntiProtonCuts) continue;
  
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
      TOF_m2_nSigma   = CalculateSigmaMassSquareTOF(Track->Pt(),TOF_m2,3,RunNumber);

    }

    double ITS_dEdx	    = -999.0;
    double ITS_dEdx_nSigma  = -999.0;
    int ITS_nCluster	    = 0;

    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);
    bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    
    if(ITSisOK){

      ITS_dEdx	      = Track->GetITSsignal();
      ITS_dEdx_nSigma = CalculateSigmadEdxITS(*Track,3,RunNumber);
      ITS_nCluster    = Track->GetITSNcls();
      if(ITS_nCluster < 0) ITS_nCluster = 0;

    }


    AntiProton_px		      = Track->Px();
    AntiProton_py		      = Track->Py();
    AntiProton_pz		      = Track->Pz();
    AntiProton_pTPC		      = Track->GetTPCmomentum();
    AntiProton_Eta		      = Track->Eta();
    AntiProton_Phi		      = Track->Phi();
    AntiProton_TPC_Chi2		      = Track->GetTPCchi2();
    AntiProton_TPC_dEdx		      = Track->GetTPCsignal();
    AntiProton_TPC_dEdx_nSigma	      = (float)fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
    AntiProton_TOF_Mass2	      = (float)TOF_m2;
    AntiProton_TOF_Mass2_nSigma	      = (float)TOF_m2_nSigma;
    AntiProton_ITS_dEdx		      = (float)ITS_dEdx;
    AntiProton_ITS_dEdx_nSigma	      = (float)ITS_dEdx_nSigma;
    AntiProton_DCAxy		      = DCAxy;
    AntiProton_DCAz		      = DCAz;
    AntiProton_TPC_nCrossedRows	      = Track->GetTPCCrossedRows();
    AntiProton_TPC_nSharedCluster     = Track->GetTPCnclsS();
    AntiProton_TPC_nFindableCluster   = Track->GetTPCNclsF();
    AntiProton_TPC_nCluster	      = Track->GetTPCNcls();
    AntiProton_ITS_nCluster	      = (unsigned short)ITS_nCluster;
    AntiProton_ID		      = track;
    AntiProton_Event_Identifier	      = EventID;
 
    fTempTree_AntiProton->Fill();
    nAntiProtonsSelected++;

  } // end of antiproton loop




  if(nAntiProtonsSelected > 1){

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

      fAntiProton_Event_Multiplicity	= Multiplicity;
      fAntiProton_Event_Centrality	= Centrality;
      fAntiProton_Event_PrimaryVertexZ  = PrimaryVertexZ;
      fAntiProton_Event_BField		= BField;

      fSaveTree_AntiProton->Fill();

    } // end of loop (copy antiprotons)


  }

  fTempTree_AntiProton->Delete();












  PostData(1,fSaveTree_Proton);
  PostData(2,fSaveTree_AntiProton);
  PostData(3,fHistoList);

} // end of UserExec















void AliAnalysisTask_pp_CreateTrees_PairsOnly::Terminate(Option_t *)
{




} // end of Terminate









// calculate the TOF beta
double AliAnalysisTask_pp_CreateTrees_PairsOnly::CalculateBetaTOF(AliAODTrack &track)
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
double AliAnalysisTask_pp_CreateTrees_PairsOnly::CalculateMassSquareTOF(AliAODTrack &track)
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








double AliAnalysisTask_pp_CreateTrees_PairsOnly::CalculateSigmaMassSquareTOF(double pT, double massSq, int ParticleSpecies, int RunNumber)
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
bool AliAnalysisTask_pp_CreateTrees_PairsOnly::CheckProtonCuts(AliAODTrack &Track, AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, bool UseOpenCuts)
{

  bool PassedParticleCuts = false;

  double Proton_pT_min, Proton_pT_max, Proton_eta_min, Proton_eta_max;
  double Proton_DCAxy_max, Proton_DCAz_max;
  double Proton_TPC_RatioRowsFindableCluster_min;
  double Proton_TPC_dEdx_nSigma_max, Proton_TPC_Chi2perCluster_max, Proton_TPC_Chi2perNDF_max;
  int Proton_TPC_nCluster_min, Proton_TPC_nCrossedRows_min, Proton_TPC_nSharedCluster_max;
  double Proton_TPC_Threshold;
  double Proton_TOF_m2_nSigma_max, Proton_TOF_m2_nSigma_max_low_pTPC;
  double Proton_ITS_dEdx_nSigma_max;
  int Proton_ITS_nCluster_min;
  bool UseTOF = true;
  bool UseITS = true;

  if(UseOpenCuts == true){

    // define open proton and antiproton track cuts
    Proton_pT_min = 0.0;
    Proton_pT_max = 4.0;
    Proton_eta_min = -0.9;
    Proton_eta_max = +0.9;
    Proton_DCAxy_max = 0.3; // cm
    Proton_DCAz_max = 0.2; // cm

    Proton_TPC_RatioRowsFindableCluster_min = 0.73;
    Proton_TPC_dEdx_nSigma_max = 4.0;
    Proton_TPC_Chi2perCluster_max = 5.0;
    Proton_TPC_Chi2perNDF_max = 5.0;
    Proton_TPC_nCluster_min = 70;
    Proton_TPC_nCrossedRows_min = 60;
    Proton_TPC_nSharedCluster_max = 2;
    Proton_TPC_Threshold = 0.8;

    Proton_TOF_m2_nSigma_max = 4.0;
    Proton_TOF_m2_nSigma_max_low_pTPC = 8.0;

    Proton_ITS_dEdx_nSigma_max = 4.0;
    Proton_ITS_nCluster_min = 1;

    UseTOF = true;
    UseITS = true;

  } // end of UseOpenCuts == true



  if(UseOpenCuts == false){

    // define closed proton and antiproton track cuts
    Proton_pT_min = 0.0;
    Proton_pT_max = 4.0;
    Proton_eta_min = -0.8;
    Proton_eta_max = +0.8;
    Proton_DCAxy_max = 0.2; // cm
    Proton_DCAz_max = 0.1; // cm

    Proton_TPC_RatioRowsFindableCluster_min = 0.83;
    Proton_TPC_dEdx_nSigma_max = 3.0;
    Proton_TPC_Chi2perCluster_max = 4.0;
    Proton_TPC_Chi2perNDF_max = 4.0;
    Proton_TPC_nCluster_min = 80;
    Proton_TPC_nCrossedRows_min = 70;
    Proton_TPC_nSharedCluster_max = 0;
    Proton_TPC_Threshold = 0.7;

    Proton_TOF_m2_nSigma_max = 3.0;
    Proton_TOF_m2_nSigma_max_low_pTPC = 7.0;

    Proton_ITS_dEdx_nSigma_max = 3.0;
    Proton_ITS_nCluster_min = 2;

    UseTOF = false;
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









































double AliAnalysisTask_pp_CreateTrees_PairsOnly::CalculateSigmadEdxITS(AliAODTrack &Track, int ParticleSpecies, int RunNumber){

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



