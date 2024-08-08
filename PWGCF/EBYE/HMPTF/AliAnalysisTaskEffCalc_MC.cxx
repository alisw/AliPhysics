
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>

#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDv0.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskEffCalc_MC.h"

//For MC event
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TList.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"




using namespace std;
using std::cout;
using std::endl;

class AliAnalysisTaskEffCalc_MC;
ClassImp(AliAnalysisTaskEffCalc_MC)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskEffCalc_MC::AliAnalysisTaskEffCalc_MC():
  AliAnalysisTaskSE(),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fMCevent(0),
  fMCstack(0),
  fInputEvent(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeList(0),
  fTreeEvent_gen(0),
  fTreeEvent_rec(0),
  fESDtrackCuts(0),
  fESDtrackCuts_primary(0),
  //fTrigger(AliVEvent::kINT7),
  fMultLow(0),
  fMultHigh(100),
  hNumberOfEvents(0),
  hNumberOfKaonPlus(0),
  hNumberOfKaonMinus(0),
  hNumberOfPionPlus(0),
  hNumberOfPionMinus(0),
  hNumberOfProtonPlus(0),
  hNumberOfProtonMinus(0),
  fTreeVariableCentrality(0),
  hist_KaonPlusWithoutPdg(0),
  hist_KaonPlusWithPdg(0),
  hist_KaonMinusWithoutPdg(0),
  hist_KaonMinusWithPdg(0),
  hist_PionPlusWithoutPdg(0),
  hist_PionPlusWithPdg(0),
  hist_PionMinusWithoutPdg(0),
  hist_PionMinusWithPdg(0),
  hist_ProtonPlusWithoutPdg(0),
  hist_ProtonPlusWithPdg(0),
  hist_ProtonMinusWithoutPdg(0),
  hist_ProtonMinusWithPdg(0),
  hist_GenKaonPlus(0),
  hist_GenKaonMinus(0),
  hist_GenProtonPlus(0),
  hist_GenProtonMinus(0),
  hist_GenPionPlus(0),
  hist_GenPionMinus(0),
  hist_GenHadronPlus(0),
  hist_GenHadronMinus(0),
  hist_RecHadronPlus(0),
  hist_RecHadronMinus(0),
  hist2D_GenProtonPlus(0),
  hist2D_GenProtonMinus(0),
  hist2D_RecProtonPlus(0),
  hist2D_RecProtonMinus(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaPionCut(0),
  fPIDnSigmaKaonCut(0),
  fPIDnSigmaProtonCut(0),
  fTPCcrossedrows(0),
  fEtaMax(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  hist_beforeCut_DCAxy(0),
  hist_beforeCut_DCAz(0),
  hist_beforeCut_eta(0),
  hist_beforeCut_chi2perTPCclstr(0),
  hist_beforeCut_chi2perITSclstr(0),
  hist_beforeCut_TPCncrossedrows(0),
  hist_afterCut_DCAxy(0),
  hist_afterCut_DCAz(0),
  hist_afterCut_eta(0),
  hist_afterCut_chi2perTPCclstr(0),
  hist_afterCut_chi2perITSclstr(0),
  hist_afterCut_TPCncrossedrows(0),
  f2Dhist_nSigmaTPC_pion(0),
  f2Dhist_nSigmaTPC_kaon(0),
  f2Dhist_nSigmaTPC_proton(0),
  f2Dhist_nSigmaTOF_pion(0),
  f2Dhist_nSigmaTOF_kaon(0),
  f2Dhist_nSigmaTOF_proton(0),
  f2Dhist_nSigmaTPCplusTOF_pion(0),
  f2Dhist_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_nSigmaTPCplusTOF_proton(0),
  fPileupCutVal(0),
  fCentralityEstimator_flag(0),
  fBayesianPID_flag(0),
  fPIDbayesPion(0),
  fPIDbayesKaon(0),
  fPIDbayesProton(0)
{
  for(int i=0; i<9; i++)
    {
      fEffProtonPlus[i] = NULL;
      fEffProtonMinus[i] = NULL;
      fEffPionPlus[i] = NULL;
      fEffPionMinus[i] = NULL;
      fEffKaonPlus[i] = NULL;
      fEffKaonMinus[i] = NULL;
    }
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskEffCalc_MC::AliAnalysisTaskEffCalc_MC(const char *name):
  AliAnalysisTaskSE(name),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fMCevent(0),
  fMCstack(0),
  fInputEvent(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeList(0),
  fTreeEvent_gen(0),
  fTreeEvent_rec(0),
  fESDtrackCuts(0),
  fESDtrackCuts_primary(0),
  //fTrigger(AliVEvent::kINT7),
  fMultLow(0),
  fMultHigh(100),
  hNumberOfEvents(0),
  hNumberOfKaonPlus(0),
  hNumberOfKaonMinus(0),
  hNumberOfPionPlus(0),
  hNumberOfPionMinus(0),
  hNumberOfProtonPlus(0),
  hNumberOfProtonMinus(0),
  fTreeVariableCentrality(0),
  hist_KaonPlusWithoutPdg(0),
  hist_KaonPlusWithPdg(0),
  hist_KaonMinusWithoutPdg(0),
  hist_KaonMinusWithPdg(0),
  hist_PionPlusWithoutPdg(0),
  hist_PionPlusWithPdg(0),
  hist_PionMinusWithoutPdg(0),
  hist_PionMinusWithPdg(0),
  hist_ProtonPlusWithoutPdg(0),
  hist_ProtonPlusWithPdg(0),
  hist_ProtonMinusWithoutPdg(0),
  hist_ProtonMinusWithPdg(0),
  hist_GenKaonPlus(0),
  hist_GenKaonMinus(0),
  hist_GenProtonPlus(0),
  hist_GenProtonMinus(0),
  hist_GenPionPlus(0),
  hist_GenPionMinus(0),
  hist_GenHadronPlus(0),
  hist_GenHadronMinus(0),
  hist_RecHadronPlus(0),
  hist_RecHadronMinus(0),
  hist2D_GenProtonPlus(0),
  hist2D_GenProtonMinus(0),
  hist2D_RecProtonPlus(0),
  hist2D_RecProtonMinus(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaPionCut(0),
  fPIDnSigmaKaonCut(0),
  fPIDnSigmaProtonCut(0),
  fTPCcrossedrows(0),
  fEtaMax(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  hist_beforeCut_DCAxy(0),
  hist_beforeCut_DCAz(0),
  hist_beforeCut_eta(0),
  hist_beforeCut_chi2perTPCclstr(0),
  hist_beforeCut_chi2perITSclstr(0),
  hist_beforeCut_TPCncrossedrows(0),
  hist_afterCut_DCAxy(0),
  hist_afterCut_DCAz(0),
  hist_afterCut_eta(0),
  hist_afterCut_chi2perTPCclstr(0),
  hist_afterCut_chi2perITSclstr(0),
  hist_afterCut_TPCncrossedrows(0),
  f2Dhist_nSigmaTPC_pion(0),
  f2Dhist_nSigmaTPC_kaon(0),
  f2Dhist_nSigmaTPC_proton(0),
  f2Dhist_nSigmaTOF_pion(0),
  f2Dhist_nSigmaTOF_kaon(0),
  f2Dhist_nSigmaTOF_proton(0),
  f2Dhist_nSigmaTPCplusTOF_pion(0),
  f2Dhist_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_nSigmaTPCplusTOF_proton(0),
  fPileupCutVal(0),
  fCentralityEstimator_flag(0),
  fBayesianPID_flag(0),
  fPIDbayesPion(0),
  fPIDbayesKaon(0),
  fPIDbayesProton(0)
{
  for(int i=0; i<9; i++)
    {
      fEffProtonPlus[i] = NULL;
      fEffProtonMinus[i] = NULL;
      fEffPionPlus[i] = NULL;
      fEffPionMinus[i] = NULL;
      fEffKaonPlus[i] = NULL;
      fEffKaonMinus[i] = NULL;
    }
 
  fUtils = new AliAnalysisUtils();
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());

}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskEffCalc_MC::~AliAnalysisTaskEffCalc_MC()  {

  if (fOutputList){
    delete fOutputList;
    fOutputList = 0x0;
  }
  if (fQAList){
    delete fQAList;
    fQAList = 0x0;
  }
  if (fTreeList){
    delete fTreeList;
    fTreeList = 0x0;
  }
  if (fTreeEvent_gen){
    delete fTreeEvent_gen;
    fTreeEvent_gen = 0x0;
  }
  if (fTreeEvent_rec){
    delete fTreeEvent_rec;
    fTreeEvent_rec = 0x0;
  }
  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }


}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskEffCalc_MC::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fTreeList   = new TList();
    fOutputList -> SetOwner(kTRUE);
    fQAList     -> SetOwner(kTRUE);
    fTreeList   -> SetOwner(kTRUE);

    OpenFile(1);
    OpenFile(2);
    OpenFile(3);
    OpenFile(4);
    OpenFile(5);

    
    //QA Plots of Event Selection
    fAODeventCuts.AddQAplotsToList(fQAList,kTRUE);
    fAODeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE, fPileupCutVal);
    
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    // hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    
    
    //distributions for calculating Kaon purity
    hist_KaonPlusWithoutPdg = new TH1F("hist_KaonPlusWithoutPdg"," hist_KaonPlusWithoutPdg",600,0,6.0);
    hist_KaonPlusWithPdg = new TH1F("hist_KaonPlusWithPdg"," hist_KaonPlusWithPdg",600,0,6.0);
    hist_KaonMinusWithoutPdg = new TH1F("hist_KaonMinusWithoutPdg"," hist_KaonMinusWithoutPdg",600,0,6.0);
    hist_KaonMinusWithPdg = new TH1F("hist_KaonMinusWithPdg"," hist_KaonMinusWithPdg",600,0,6.0);
    fOutputList->Add(hist_KaonPlusWithoutPdg);
    fOutputList->Add(hist_KaonPlusWithPdg);
    fOutputList->Add(hist_KaonMinusWithoutPdg);
    fOutputList->Add(hist_KaonMinusWithPdg);

    //distributions for calculating Pion purity
    hist_PionPlusWithoutPdg = new TH1F("hist_PionPlusWithoutPdg"," hist_PionPlusWithoutPdg",600,0,6.0);
    hist_PionPlusWithPdg = new TH1F("hist_PionPlusWithPdg"," hist_PionPlusWithPdg",600,0,6.0);
    hist_PionMinusWithoutPdg = new TH1F("hist_PionMinusWithoutPdg"," hist_PionMinusWithoutPdg",600,0,6.0);
    hist_PionMinusWithPdg = new TH1F("hist_PionMinusWithPdg"," hist_PionMinusWithPdg",600,0,6.0);
    fOutputList->Add(hist_PionPlusWithoutPdg);
    fOutputList->Add(hist_PionPlusWithPdg);
    fOutputList->Add(hist_PionMinusWithoutPdg);
    fOutputList->Add(hist_PionMinusWithPdg);

    //distributions for calculating Proton purity
    hist_ProtonPlusWithoutPdg = new TH1F("hist_ProtonPlusWithoutPdg"," hist_ProtonPlusWithoutPdg",600,0,6.0);
    hist_ProtonPlusWithPdg = new TH1F("hist_ProtonPlusWithPdg"," hist_ProtonPlusWithPdg",600,0,6.0);
    hist_ProtonMinusWithoutPdg = new TH1F("hist_ProtonMinusWithoutPdg"," hist_ProtonMinusWithoutPdg",600,0,6.0);
    hist_ProtonMinusWithPdg = new TH1F("hist_ProtonMinusWithPdg"," hist_ProtonMinusWithPdg",600,0,6.0);
    fOutputList->Add(hist_ProtonPlusWithoutPdg);
    fOutputList->Add(hist_ProtonPlusWithPdg);
    fOutputList->Add(hist_ProtonMinusWithoutPdg);
    fOutputList->Add(hist_ProtonMinusWithPdg);

    //Generated distributions: for efficiency
    hist_GenKaonPlus = new TH1F("hist_GenKaonPlus"," hist_GenKaonPlus",600,0,6.0);
    hist_GenKaonMinus = new TH1F("hist_GenKaonMinus"," hist_GenKaonMinus",600,0,6.0);
    hist_GenPionPlus = new TH1F("hist_GenPionPlus"," hist_GenPionPlus",600,0,6.0);
    hist_GenPionMinus = new TH1F("hist_GenPionMinus"," hist_GenPionMinus",600,0,6.0);
    hist_GenProtonPlus = new TH1F("hist_GenProtonPlus"," hist_GenProtonPlus",600,0,6.0);
    hist_GenProtonMinus = new TH1F("hist_GenProtonMinus"," hist_GenProtonMinus",600,0,6.0);
    hist_GenHadronPlus = new TH1F("hist_GenHadronPlus"," hist_GenHadronPlus",600,0,6.0);
    hist_GenHadronMinus = new TH1F("hist_GenHadronMinus"," hist_GenHadronMinus",600,0,6.0);
    fOutputList->Add(hist_GenKaonPlus);
    fOutputList->Add(hist_GenKaonMinus);
    fOutputList->Add(hist_GenPionPlus);
    fOutputList->Add(hist_GenPionMinus);
    fOutputList->Add(hist_GenProtonPlus);
    fOutputList->Add(hist_GenProtonMinus);
    fOutputList->Add(hist_GenHadronPlus);
    fOutputList->Add(hist_GenHadronMinus);

    hist_RecHadronPlus = new TH1F("hist_RecHadronPlus"," hist_RecHadronPlus",600,0,6.0);
    hist_RecHadronMinus = new TH1F("hist_RecHadronMinus"," hist_RecHadronMinus",600,0,6.0);
    fOutputList->Add(hist_RecHadronPlus);
    fOutputList->Add(hist_RecHadronMinus);


    hist2D_GenProtonPlus = new TH2D("hist2D_GenProtonPlus", "hist2D_GenProtonPlus",100,0,100,600,0,6.0);
    hist2D_GenProtonMinus  = new TH2D("hist2D_GenProtonMinus", "hist2D_GenProtonMinus",100,0,100,600,0,6.0);
    hist2D_RecProtonPlus = new TH2D("hist2D_RecProtonPlus", "hist2D_RecProtonPlus",100,0,100,600,0,6.0);
    hist2D_RecProtonMinus  = new TH2D("hist2D_RecProtonMinus", "hist2D_RecProtonMinus",100,0,100,600,0,6.0);
    fOutputList->Add(hist2D_GenProtonPlus);
    fOutputList->Add(hist2D_GenProtonMinus);
    fOutputList->Add(hist2D_RecProtonPlus);
    fOutputList->Add(hist2D_RecProtonMinus);


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Histograms for track variables

    //before track cut
    hist_beforeCut_DCAxy = new TH1D("hist_beforeCut_DCAxy","hist_beforeCut_DCAxy", 1000, -5, +5);
    hist_beforeCut_DCAz = new TH1D("hist_beforeCut_DCAz","hist_beforeCut_DCAz", 1000, -5, +5);
    hist_beforeCut_eta = new TH1D ("hist_beforeCut_eta","hist_beforeCut_eta", 20, -1, +1);
    hist_beforeCut_chi2perTPCclstr = new TH1D ("hist_beforeCut_chi2perTPCclstr", "hist_beforeCut_chi2perTPCclstr",100, 0, 5);
    hist_beforeCut_chi2perITSclstr = new TH1D ("hist_beforeCut_chi2perITSclstr", "hist_beforeCut_chi2perITSclstr",100, 0, 50);
    hist_beforeCut_TPCncrossedrows = new TH1D ("hist_beforeCut_TPCncrossedrows", "hist_beforeCut_TPCncrossedrows",200, 0, 200);
    fOutputList->Add(hist_beforeCut_DCAxy);
    fOutputList->Add(hist_beforeCut_DCAz);
    fOutputList->Add(hist_beforeCut_eta);
    fOutputList->Add(hist_beforeCut_chi2perTPCclstr);
    fOutputList->Add(hist_beforeCut_chi2perITSclstr);
    fOutputList->Add(hist_beforeCut_TPCncrossedrows);

    //after track cut
    hist_afterCut_DCAxy = new TH1D("hist_afterCut_DCAxy","hist_afterCut_DCAxy", 1000, -5, +5);
    hist_afterCut_DCAz = new TH1D("hist_afterCut_DCAz","hist_afterCut_DCAz", 1000, -5, +5);
    hist_afterCut_eta = new TH1D ("hist_afterCut_eta","hist_afterCut_eta", 20, -1, +1);
    hist_afterCut_chi2perTPCclstr = new TH1D ("hist_afterCut_chi2perTPCclstr", "hist_afterCut_chi2perTPCclstr",100, 0, 5);
    hist_afterCut_chi2perITSclstr = new TH1D ("hist_afterCut_chi2perITSclstr", "hist_afterCut_chi2perITSclstr",100, 0, 50);
    hist_afterCut_TPCncrossedrows = new TH1D ("hist_afterCut_TPCncrossedrows", "hist_afterCut_TPCncrossedrows",200, 0, 200);
    fOutputList->Add(hist_afterCut_DCAxy);
    fOutputList->Add(hist_afterCut_DCAz);
    fOutputList->Add(hist_afterCut_eta);
    fOutputList->Add(hist_afterCut_chi2perTPCclstr);
    fOutputList->Add(hist_afterCut_chi2perITSclstr);
    fOutputList->Add(hist_afterCut_TPCncrossedrows);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //PID nSigma Histograms 2d

    f2Dhist_nSigmaTPC_pion = new TH2D("f2Dhist_nSigmaTPC_pion", "f2Dhist_nSigmaTPC_pion", 2000, -10, +10, 500, 0, 5);
    f2Dhist_nSigmaTPC_kaon = new TH2D("f2Dhist_nSigmaTPC_kaon", "f2Dhist_nSigmaTPC_kaon", 2000, -10, +10, 500, 0, 5);
    f2Dhist_nSigmaTPC_proton = new TH2D("f2Dhist_nSigmaTPC_proton", "f2Dhist_nSigmaTPC_proton", 2000, -10, +10, 500, 0, 5);

    f2Dhist_nSigmaTOF_pion = new TH2D("f2Dhist_nSigmaTOF_pion", "f2Dhist_nSigmaTOF_pion", 2000, -10, +10, 500, 0, 5);
    f2Dhist_nSigmaTOF_kaon = new TH2D("f2Dhist_nSigmaTOF_kaon", "f2Dhist_nSigmaTOF_kaon", 2000, -10, +10, 500, 0, 5);
    f2Dhist_nSigmaTOF_proton = new TH2D("f2Dhist_nSigmaTOF_proton", "f2Dhist_nSigmaTOF_proton", 2000, -10, +10, 500, 0, 5);

     f2Dhist_nSigmaTPCplusTOF_pion = new TH2D("f2Dhist_nSigmaTPCplusTOF_pion", "f2Dhist_nSigmaTPCplusTOF_pion", 2000, -10, +10, 500, 0, 5);
    f2Dhist_nSigmaTPCplusTOF_kaon = new TH2D("f2Dhist_nSigmaTPCplusTOF_kaon", "f2Dhist_nSigmaTPCplusTOF_kaon", 2000, -10, +10, 500, 0, 5);
    f2Dhist_nSigmaTPCplusTOF_proton = new TH2D("f2Dhist_nSigmaTPCplusTOF_proton", "f2Dhist_nSigmaTPCplusTOF_proton", 2000, -10, +10, 500, 0, 5);
    fOutputList->Add(f2Dhist_nSigmaTPC_pion);
    fOutputList->Add(f2Dhist_nSigmaTPC_kaon);
    fOutputList->Add(f2Dhist_nSigmaTPC_proton);
    fOutputList->Add(f2Dhist_nSigmaTOF_pion);
    fOutputList->Add(f2Dhist_nSigmaTOF_kaon);
    fOutputList->Add(f2Dhist_nSigmaTOF_proton);
    fOutputList->Add(f2Dhist_nSigmaTPCplusTOF_pion);
    fOutputList->Add(f2Dhist_nSigmaTPCplusTOF_kaon);
    fOutputList->Add(f2Dhist_nSigmaTPCplusTOF_proton);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //TTree object to store variables
    fTreeEvent_gen = new TTree("fTreeEvent_gen","Event Tree");
    fTreeEvent_gen->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    //Generated
    
    //TTree object to store variables
    fTreeEvent_rec = new TTree("fTreeEvent_rec","Event Tree");
    fTreeEvent_rec->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    //Reconsctructed
    
   
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeList);
    PostData(4, fTreeEvent_gen);
    PostData(5, fTreeEvent_rec);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskEffCalc_MC::UserExec(Option_t *)  {
  
    //Get Input Event
    if ( !GetEvent ()) return;
    cout<<"*********************** Found AOD event !!! ******************************"<<endl;

    
    //Get multiplicity percentile
    Float_t lV0M;
    AliMultSelection *MultSelection = (AliMultSelection*) fAODevent -> FindListObject("MultSelection");
    if( !MultSelection)
      {
	return;
      }
    else
      {
	if (fCentralityEstimator_flag == 0)
	  lV0M = MultSelection->GetMultiplicityPercentile("V0M");
	else if (fCentralityEstimator_flag == 1)
	  lV0M = MultSelection->GetMultiplicityPercentile("CL0");
	else if (fCentralityEstimator_flag == 2)
	  lV0M = MultSelection->GetMultiplicityPercentile("CL1");
	else if (fCentralityEstimator_flag == 3)
	  lV0M = MultSelection->GetMultiplicityPercentile("CL2");

	//cout<<"V0M: "<<lV0M<<endl;
      }

    fMCevent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCevent) {
      Printf("ERROR: Could not retrieve MC event \n");
      cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
      PostData(4, fTreeEvent_gen);
      PostData(5, fTreeEvent_rec);
      return;
    }
  
    

    //"fTreeEvent" Tree Variable
    fTreeVariableCentrality = lV0M;

   
    //Loop on generated MC tracks
    Int_t noGenMCtracks = fMCevent->GetNumberOfTracks();
    cout<<"No of generated MC tracks: "<<noGenMCtracks<<endl;
    for(Int_t itr_mcgen=0; itr_mcgen < noGenMCtracks; itr_mcgen++)
      {
	AliAODMCParticle *mcGenTrack = (AliAODMCParticle*) fMCevent->GetTrack(itr_mcgen);
	if (!mcGenTrack)
	  {
	    cout<<"Could not find track in MC generated loop !!!"<<endl;
	    continue;
	  }
	//cout<<"Found generated MC track !!"<<endl;
	
	//Out of bunch pileup event removal
	if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(itr_mcgen,fMCevent))
	  {
	    //cout<<"Track belongs to out of bunch pileup events. Removed !!!!"<<endl;
	    continue;
	  }
  	
	Double_t TrackGen_eta = mcGenTrack->Eta();
	Double_t TrackGen_pt = mcGenTrack->Pt();
	Double_t TrackGen_charge = mcGenTrack->Charge();
	Double_t TrackGen_PID = mcGenTrack->GetPdgCode();
	Double_t Track_MotherLabel = mcGenTrack->GetMother();

	//cout<<"Found generated track after imposing all cuts !!!-------------->"<<endl;

	
	if(!mcGenTrack->IsPhysicalPrimary())  // removed IsPhysicalPrimary() condition 
	  continue;
	

	//if(Track_MotherLabel > 0) continue;   /* removed MotherLabel condition */

	if (!TrackGen_PID) continue;
	if (TrackGen_pt < 0.2) continue;
	if (TrackGen_pt > 6.0) continue;
	if (TMath::Abs(TrackGen_eta) > 0.8) continue;

	
	if(TrackGen_charge > 0)
	  hist_GenHadronPlus->Fill(TrackGen_pt);
	if(TrackGen_charge < 0)
	  hist_GenHadronMinus->Fill(TrackGen_pt);

	if (TMath::Abs(TrackGen_PID) == 321)
	  {
	    if(TrackGen_charge > 0)
	      hist_GenKaonPlus->Fill(TrackGen_pt);
	    if(TrackGen_charge < 0)
	      hist_GenKaonMinus->Fill(TrackGen_pt);
	  }

	if (TMath::Abs(TrackGen_PID) == 2212)
	  {
	    if(TrackGen_charge > 0)
	      {
		hist_GenProtonPlus->Fill(TrackGen_pt);
		hist2D_GenProtonPlus->Fill(lV0M, TrackGen_pt);
	      }
	    if(TrackGen_charge < 0)
	      {
		hist_GenProtonMinus->Fill(TrackGen_pt);
		hist2D_GenProtonMinus->Fill(lV0M, TrackGen_pt);
	      }
	  }

	if (TMath::Abs(TrackGen_PID) == 211)
	  {
	    if(TrackGen_charge > 0)
	      hist_GenPionPlus->Fill(TrackGen_pt);
	    if(TrackGen_charge < 0)
	      hist_GenPionMinus->Fill(TrackGen_pt);
	  }

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Calculation of variables for analysis using rec. tracks
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//to be filled
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      }
    //end generated track loop

    /*
    Int_t ptBinNo = 0;
    Double_t BinCont = 0;
    Double_t EffWgt = 0;

    if(fListTRKCorr) GetMCEffCorrectionHist();

    Int_t centrality_bin = 0;
    if(lV0M > 0.0 && lV0M < 5.0)
      centrality_bin = 0;
    else if (lV0M > 5.0 && lV0M < 10.0)
      centrality_bin = 1;
    else if (lV0M > 10.0 && lV0M < 20.0)
      centrality_bin = 2;
    else if (lV0M > 20.0 && lV0M < 30.0)
      centrality_bin = 3;
    else if (lV0M > 30.0 && lV0M < 40.0)
      centrality_bin = 4;
    else if (lV0M > 40.0 && lV0M < 50.0)
      centrality_bin = 5;
    else if (lV0M > 50.0 && lV0M < 60.0)
      centrality_bin = 6;
    else if (lV0M > 60.0 && lV0M < 70.0)
      centrality_bin = 7;
    else if (lV0M > 70.0 && lV0M < 90.0)
      centrality_bin = 8;
    */


    //-------------------------------------------------------------------------------------------
    //Reconstructed variables Initialization
    //-------------------------------------------------------------------------------------------
   

    //Loop on reconstructed MC tracks
    Int_t noRecMCtracks = fAODevent->GetNumberOfTracks();
    cout<<"No of reconstructed MC tracks: "<<noRecMCtracks<<endl;
    for(Int_t itr_mcrec=0; itr_mcrec < noRecMCtracks; itr_mcrec++)
      {
	
	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr_mcrec);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;

	Int_t aodtrk_label = TMath::Abs(aodtrack->GetLabel());
	AliAODMCParticle *mcRecTrack = (AliAODMCParticle*) fMCevent->GetTrack(aodtrk_label);
	if (!mcRecTrack)
	  {
	    cout<<"Could not find MC generted track for the AODtrack Label !!!"<<endl;
	    continue;
	  }
	//cout<<"Found reconstructed MC track (matched with generated) !!"<<endl;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Track parameters
	Double_t Track_pt = mcRecTrack->Pt();
	Double_t Track_charge = mcRecTrack->Charge();
	Double_t Track_PID = mcRecTrack->GetPdgCode();
	Double_t Track_eta = mcRecTrack->Eta();
	Double_t Track_MotherLabel = mcRecTrack->GetMother();
	Double_t trkITSchi2 = aodtrack->GetITSchi2();
	Int_t trkITSNcls = aodtrack->GetITSNcls();
	Double_t trkITSchi2perNcls = trkITSchi2/trkITSNcls;
	Double_t trkTPCchi2perNcls = aodtrack->GetTPCchi2perCluster();
	Double_t trkTPCcrossedrows = aodtrack->GetTPCCrossedRows();

	//track DCAxy and DCAz
	Float_t trkDCAxy;
	Float_t trkDCAz;
	aodtrack->GetImpactParameters(trkDCAxy, trkDCAz);
	
	
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Histograms filled befor applying track cut
	hist_beforeCut_DCAxy->Fill(trkDCAxy);
	hist_beforeCut_DCAz->Fill(trkDCAz);
	hist_beforeCut_eta->Fill(Track_eta);
	hist_beforeCut_chi2perTPCclstr->Fill(trkTPCchi2perNcls);
	hist_beforeCut_chi2perITSclstr->Fill(trkITSchi2perNcls);
	hist_beforeCut_TPCncrossedrows->Fill(trkTPCcrossedrows);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Track quality cuts as in data
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Filterbit
	if(!aodtrack->TestFilterBit(fFBNo))  continue; /*default: 96, systematics: 768*/

	//cuts on TPCchi2perClstr and ITSchi2perClstr and TPCnCrossedRows
	if (trkTPCcrossedrows < fTPCcrossedrows) continue;
	if (trkTPCchi2perNcls > fChi2TPC) continue;
	if (trkITSchi2perNcls > fChi2ITS) continue;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	

	if(!mcRecTrack->IsPhysicalPrimary()) continue;

	//if(Track_MotherLabel > 0) continue;

	if (!Track_PID) continue;

	//cout<<"Found a reconstructed track!!"<<endl;
	/*
	//########################################################################################
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//PID nSigma histograms
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Double_t fTPCnSigma_Pion = 0.0;
	Double_t fTPCnSigma_Proton = 0.0;
	Double_t fTPCnSigma_Kaon = 0.0;
	Double_t fTOFnSigma_Pion = 0.0;
	Double_t fTOFnSigma_Proton = 0.0;
	Double_t fTOFnSigma_Kaon = 0.0;

	//TPC nsigma
	fTPCnSigma_Pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
	fTPCnSigma_Proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
	fTPCnSigma_Kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
	//TOF nsigma
	fTOFnSigma_Pion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
	fTOFnSigma_Proton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
	fTOFnSigma_Kaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

	Double_t fTPCplusTOFnSigma_Pion = sqrt(TMath::Power(fTPCnSigma_Pion, 2.0) + TMath::Power(fTOFnSigma_Pion, 2.0));
	Double_t fTPCplusTOFnSigma_Kaon = sqrt(TMath::Power(fTPCnSigma_Kaon, 2.0) + TMath::Power(fTOFnSigma_Kaon, 2.0));
	Double_t fTPCplusTOFnSigma_Proton = sqrt(TMath::Power(fTPCnSigma_Proton, 2.0) + TMath::Power(fTOFnSigma_Proton, 2.0));

	//TPC
	f2Dhist_nSigmaTPC_pion->Fill(fTPCnSigma_Pion, Track_pt);
	f2Dhist_nSigmaTPC_kaon->Fill(fTPCnSigma_Kaon, Track_pt);
	f2Dhist_nSigmaTPC_proton->Fill(fTPCnSigma_Proton, Track_pt);
	//TOF
	f2Dhist_nSigmaTOF_pion->Fill(fTOFnSigma_Pion, Track_pt);
	f2Dhist_nSigmaTOF_kaon->Fill(fTOFnSigma_Kaon, Track_pt);
	f2Dhist_nSigmaTOF_proton->Fill(fTOFnSigma_Proton, Track_pt);
	//TPC+TOF
	f2Dhist_nSigmaTPCplusTOF_pion->Fill(fTPCplusTOFnSigma_Pion, Track_pt);
	f2Dhist_nSigmaTPCplusTOF_kaon->Fill(fTPCplusTOFnSigma_Kaon, Track_pt);
	f2Dhist_nSigmaTPCplusTOF_proton->Fill(fTPCplusTOFnSigma_Proton, Track_pt);

	//########################################################################################
	*/
	
	//Kinematic acceptance cut for analysis
	if (Track_pt < 0.2) continue;
	if (Track_pt > 6.0) continue;
	if (TMath::Abs(Track_eta) > 0.8) continue;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Histograms filled after applying track cut
	hist_afterCut_DCAxy->Fill(trkDCAxy);
	hist_afterCut_DCAz->Fill(trkDCAz);
	hist_afterCut_eta->Fill(Track_eta);
	hist_afterCut_chi2perTPCclstr->Fill(trkTPCchi2perNcls);
	hist_afterCut_chi2perITSclstr->Fill(trkITSchi2perNcls);
	hist_afterCut_TPCncrossedrows->Fill(trkTPCcrossedrows);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if(Track_charge > 0)
	  hist_RecHadronPlus->Fill(Track_pt);
	if(Track_charge < 0)
	  hist_RecHadronMinus->Fill(Track_pt);

	//PID selection
	Bool_t IsPion = kFALSE;
	Bool_t IsKaon = kFALSE;
	Bool_t IsProton = kFALSE;

	if(fBayesianPID_flag == 1)
	  {
	    Int_t trackPIDbasedId = IdentifyTrackBayesian(track);
	    if (trackPIDbasedId == 1)
	      {
		IsPion = kTRUE;
	      }
	    else if (trackPIDbasedId == 2)
	      {
		IsKaon = kTRUE; 
	      }
	    else if(trackPIDbasedId == 3)
	      {
		IsProton = kTRUE;
	      }
	    else
	      {
		IsPion = kFALSE; IsKaon = kFALSE; IsProton = kFALSE;
	      }
	  }
	else //Traditional method of PID by nSigma
	  {
	    IsKaon = KaonSelector(track, fPIDnSigmaKaonCut);
	    IsPion = PionSelector(track, fPIDnSigmaPionCut);
	    IsProton = ProtonSelector(track, fPIDnSigmaProtonCut);
	  }

	
	if (!IsKaon && !IsPion && !IsProton) continue;

	/*
	Int_t flag = 0;
	if(IsKaon) flag+=1;
	if(IsPion) flag+=1;
	if(IsProton) flag+=1;
	//cout<<"Particle identified as more than on PID: flag= "<<flag<<endl;
	if(flag>1) continue;
	*/

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Filling histograms to calculate purity
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	if (Track_charge > 0 && IsKaon)   //K+
	  {
	    //Histograms for purity
	    
	    hist_KaonPlusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 321)
	      hist_KaonPlusWithPdg->Fill(Track_pt);
	  }
	  

	if (Track_charge < 0 && IsKaon)   //K-
	  {
	    //Histograms for purity
	    
	    hist_KaonMinusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 321)
	      hist_KaonMinusWithPdg->Fill(Track_pt);
	  }

	if (Track_charge > 0 && IsProton)   //proton
	  {
	    //Histograms for purity
	    
	    hist_ProtonPlusWithoutPdg->Fill(Track_pt);
	    hist2D_RecProtonPlus->Fill(lV0M, Track_pt);
	    if (TMath::Abs(Track_PID) == 2212)
	      hist_ProtonPlusWithPdg->Fill(Track_pt);
	  }
	  

	if (Track_charge < 0 && IsProton)   //anti-proton
	  {
	    //Histograms for purity
	    
	    hist_ProtonMinusWithoutPdg->Fill(Track_pt);
	    hist2D_RecProtonMinus->Fill(lV0M, Track_pt);
	    if (TMath::Abs(Track_PID) == 2212)
	      hist_ProtonMinusWithPdg->Fill(Track_pt);
	  }

	if (Track_charge > 0 && IsPion)   //pi+
	  {
	    //Histograms for purity
	    
	    hist_PionPlusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 211)
	      hist_PionPlusWithPdg->Fill(Track_pt);
	  }
	  

	if (Track_charge < 0 && IsPion)   //pi-
	  {
	    //Histograms for purity
	    
	    hist_PionMinusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 211)
	      hist_PionMinusWithPdg->Fill(Track_pt);
	  }


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Calculation of variables for analysis using rec. tracks
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//to be filled
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	
      }
    //end reconstructed track loop

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Filling tree variables
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    fTreeVariableCentrality=lV0M;
    fTreeEvent_gen->Fill();
    fTreeEvent_rec->Fill();
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeList);
    PostData(4, fTreeEvent_gen);
    PostData(5, fTreeEvent_rec);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::GetEvent ()  //event cuts copied from my code written earlier 

{
 
  //Get Input Event
  fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
  if (!fAODevent)
    {
      AliWarning("ERROR: lAODevent not available from InputEvent(), trying with AODEvent() \n");

      fAODevent  = AODEvent();
      if(!fAODevent)
	{
	  AliWarning("ERROR: lAODevent not available from AODEvent() Aborting event!");
	  PostData(1, fOutputList);
	  PostData(2, fQAList);
	  PostData(3, fTreeList);
	  PostData(4, fTreeEvent_gen);
	  PostData(5, fTreeEvent_rec);
	  return kFALSE;
	}
    }

  fInputEvent = dynamic_cast <AliVEvent*>(InputEvent());
  if(!fInputEvent)
    {
      AliWarning("ERROR: fInputEvent (AliVEvent) not available \n");
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
      PostData(4, fTreeEvent_gen);
      PostData(5, fTreeEvent_rec);
      return kFALSE;
    }
  
  hNumberOfEvents -> Fill(0.5);

    
  //Standard Event Cuts
  if (!fAODeventCuts.AcceptEvent(fInputEvent)) {
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeList);
    PostData(4, fTreeEvent_gen);
    PostData(5, fTreeEvent_rec);
    return kFALSE;
  }
  hNumberOfEvents -> Fill(1.5);
  
  
  //Reject Events with Incomplete DAQ
  if (fAODevent->IsIncompleteDAQ())
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
      PostData(4, fTreeEvent_gen);
      PostData(5, fTreeEvent_rec);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(2.5);
        
  //V0 Timing Decision
  AliVVZERO *vzeroData = fInputEvent->GetVZEROData();
  if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision()))
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
      PostData(4, fTreeEvent_gen);
      PostData(5, fTreeEvent_rec);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(3.5);
  
        
  //Pileup Rejection
  Int_t nClustersLayer0 = fInputEvent->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = fInputEvent->GetNumberOfITSClusters(1);
  Int_t nTracklets      = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
  if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0)
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
      PostData(4, fTreeEvent_gen);
      PostData(5, fTreeEvent_rec);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(4.5);

    
  //Primary Vertex Tracks
  AliAODVertex *vertex_tracks = (AliAODVertex*) fAODevent->GetPrimaryVertexTracks();
  if (!vertex_tracks)
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
      PostData(4, fTreeEvent_gen);
      PostData(5, fTreeEvent_rec);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(5.5);
  
        
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(6.5);

    
        
    //Primary Vertex SPD
    AliAODVertex *vertex_SPD = (AliAODVertex*) fAODevent->GetPrimaryVertexSPD();
    if (!vertex_SPD)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(7.5);
        
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(8.5);

    
    //SPD Pile-up in Mult Bins
    if (fAODevent->IsPileupFromSPDInMultBins())
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(9.5);

    
    
    ////tigger/////////////                                                                                                                     
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
    if ( !isSelected)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(10.5);
    
        
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(11.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -1.0*fVertexZMax || vertex_tracks->GetZ() > +1.0*fVertexZMax)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(12.5);
               
    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      } 
    hNumberOfEvents -> Fill(13.5);
                
    //Selection of Multiplicity Range
    Double_t mult_percentile;
    if (fCentralityEstimator_flag == 0)
      mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    else if (fCentralityEstimator_flag == 1)
      mult_percentile = multiplicitySelection->GetMultiplicityPercentile("CL0");
    else if (fCentralityEstimator_flag == 2)
      mult_percentile = multiplicitySelection->GetMultiplicityPercentile("CL1");
    else if (fCentralityEstimator_flag == 3)
      mult_percentile = multiplicitySelection->GetMultiplicityPercentile("CL2");
    if (mult_percentile < 0.0 || mult_percentile > 90.0)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	PostData(4, fTreeEvent_gen);
	PostData(5, fTreeEvent_rec);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(14.5);
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();
    fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask
    
    return kTRUE;
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::PassedTrackQualityCuts (AliAODTrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);

    //select global tracks: ITS+TPC loose DCA
    // if(!track->TestFilterMask(BIT(4)))
    //   return passedTrkSelection;
    
    //Track Selection Cuts
    Int_t nTPCcluster = track->GetTPCNcls();  //TPC cluster cut
    if (nTPCcluster < 70)
      return passedTrkSelection;

    Double_t chi2TPCperClstr = track->GetTPCchi2perCluster();
    if(chi2TPCperClstr > 4)
      return passedTrkSelection;

    
    if(!(AliAODTrack::kTPCrefit)) //TPC refit
      return passedTrkSelection;

    if(track->GetKinkIndex(0) > 0) //No kink daughters
      return passedTrkSelection;

    if(!(AliAODTrack::kITSrefit)) //ITS refit
      return passedTrkSelection;
    
    
    //Kinematic cuts
    Double_t eta = track->Eta(); //Eta cut
    if(TMath::Abs(eta) > 0.8)
      return passedTrkSelection;

    Double_t pt = track->Pt(); //Pt cut
    if(pt < 0.15)
      return passedTrkSelection;

    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
 //_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::KaonSelector(AliVTrack *track,  Double_t nSigmaKaonCut)  {
 
  Double_t p[3];
  track->PxPyPz(p);

  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fTPCnSigmaPion = 0.0;
  Double_t fTPCnSigmaProton = 0.0;
  Double_t fTPCnSigmaKaon = 0.0;
  Double_t fTOFnSigmaPion = 0.0;
  Double_t fTOFnSigmaProton = 0.0;
  Double_t fTOFnSigmaKaon = 0.0;

  //TPC nsigma
  fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //TOF nsigma
  fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
  Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
  Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

 
  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.5 )
    {
      /*
      //rejection
      Int_t flag = 0;
      if (TMath::Abs(fTPCnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCnSigmaKaon > fTPCnSigmaProton) return kFALSE;
      if (fTPCnSigmaKaon > fTPCnSigmaPion) return kFALSE;
      */
      
      //acception
      if(TMath::Abs(fTPCnSigmaKaon) < nSigmaKaonCut)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.5 )
    {
      //rejection
      Int_t flag = 0;
      if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCplusTOFnSigmaKaon > fTPCplusTOFnSigmaProton) return kFALSE;
      if (fTPCplusTOFnSigmaKaon > fTPCplusTOFnSigmaPion) return kFALSE;

      //acception
      if (fTPCplusTOFnSigmaKaon < nSigmaKaonCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
  
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::PionSelector(AliVTrack *track,  Double_t nSigmaPionCut)  {
  
  Double_t p[3];
  track->PxPyPz(p);

  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fTPCnSigmaPion = 0.0;
  Double_t fTPCnSigmaProton = 0.0;
  Double_t fTPCnSigmaKaon = 0.0;
  Double_t fTOFnSigmaPion = 0.0;
  Double_t fTOFnSigmaProton = 0.0;
  Double_t fTOFnSigmaKaon = 0.0;

  //TPC nsigma
  fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //TOF nsigma
  fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
  Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
  Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

  


  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.5 )
    {
      /*
      //rejection
      
      Int_t flag = 0;
      if (TMath::Abs(fTPCnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCnSigmaPion > fTPCnSigmaProton) return kFALSE;
      if (fTPCnSigmaPion > fTPCnSigmaKaon) return kFALSE;
      */
      
      //acception
      
      if(TMath::Abs(fTPCnSigmaPion) < nSigmaPionCut)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.5 )
    {
      //rejection
      
      Int_t flag = 0;
      if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCplusTOFnSigmaPion > fTPCplusTOFnSigmaProton) return kFALSE;
      if (fTPCplusTOFnSigmaPion > fTPCplusTOFnSigmaKaon) return kFALSE;

      //acception
      
      if (fTPCplusTOFnSigmaPion < nSigmaPionCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::ProtonSelector(AliVTrack *track,  Double_t nSigmaProtonCut)  {
  
  Double_t p[3];
  track->PxPyPz(p);

  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fTPCnSigmaPion = 0.0;
  Double_t fTPCnSigmaProton = 0.0;
  Double_t fTPCnSigmaKaon = 0.0;
  Double_t fTOFnSigmaPion = 0.0;
  Double_t fTOFnSigmaProton = 0.0;
  Double_t fTOFnSigmaKaon = 0.0;

  //TPC nsigma
  fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //TOF nsigma
  fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
  Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
  Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

  
  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.6 )
    {
      //rejection

      /*
      Int_t flag = 0;
      if (TMath::Abs(fTPCnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCnSigmaProton > fTPCnSigmaPion) return kFALSE;
      if (fTPCnSigmaProton > fTPCnSigmaKaon) return kFALSE;
      */
      
      //acception

      if(TMath::Abs(fTPCnSigmaProton) < nSigmaProtonCut)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.6 )
    {
      //rejection

      
      Int_t flag = 0;
      if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCplusTOFnSigmaProton > fTPCplusTOFnSigmaPion) return kFALSE;
      if (fTPCplusTOFnSigmaProton > fTPCplusTOFnSigmaKaon) return kFALSE;
      

      //acception
      
      if (fTPCplusTOFnSigmaProton < nSigmaProtonCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
}


//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
    //Initialization
    Bool_t passedPIDSelection=(kFALSE);
    
    //TPC Particle Identification
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,type);
    if (nsigmaTPC < -4.0) return passedPIDSelection;
    if (nsigmaTPC > +4.0) return passedPIDSelection;

    passedPIDSelection = kTRUE;
    return passedPIDSelection;
}
 //_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kTRUE);
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  return passedTrackPileupCut;
}
 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskEffCalc_MC::GetMCEffCorrectionHist()
{
  if(fListTRKCorr)
    {
      fHistMCEffKaonPlus = (TH1D*) fListTRKCorr->FindObject("histKaonPlusEff");
      fHistMCEffKaonMinus = (TH1D*) fListTRKCorr->FindObject("histKaonMinusEff");

      fHistMCEffPionPlus = (TH1D*) fListTRKCorr->FindObject("histPionPlusEff");
      fHistMCEffPionMinus = (TH1D*) fListTRKCorr->FindObject("histPionMinusEff");

      fHistMCEffProtonPlus = (TH1D*) fListTRKCorr->FindObject("histProtonPlusEff");
      fHistMCEffProtonMinus = (TH1D*) fListTRKCorr->FindObject("histProtonMinusEff");

      for(int i=0; i<9; i++)
	{
	  fEffProtonPlus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffProtonPlus%d",i));
	  fEffProtonMinus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffProtonMinus%d",i));
	  fEffPionPlus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffPionPlus%d",i));
	  fEffPionMinus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffPionMinus%d",i));
	  fEffKaonPlus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffKaonPlus%d",i));
	  fEffKaonMinus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffKaonMinus%d",i));

	}
    }

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::HasTrackPIDTPC(AliVTrack *track)
{
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  return (pidStatusTPC == AliPIDResponse::kDetPidOk);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskEffCalc_MC::HasTrackPIDTOF(AliVTrack *track) 
{
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  return ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskEffCalc_MC::IdentifyTrackBayesian(AliVTrack *track) // identify Pi, Ka, Pr based on BayesianPID
{
  // checking detector statuses
  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  if(!bIsTPCok) { return -1; }
  

  Double_t l_Probs[AliPID::kSPECIES];
  Double_t l_MaxProb[] = {fPIDbayesPion,fPIDbayesKaon,fPIDbayesProton};
  
  UInt_t flag=fPIDCombined->ComputeProbabilities(track, fPIDResponse, l_Probs);
  Bool_t l_TOFUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, l_Probs) & AliPIDResponse::kDetTOF;
  Int_t pidInd = 0;
  for(Int_t i(0); i < AliPID::kSPECIES; i++)
    pidInd=(l_Probs[i]>l_Probs[pidInd])?i:pidInd;
  Int_t retInd = pidInd-AliPID::kPion+1; //realigning
  if(retInd<1 || retInd>3) return -1;
  if(l_Probs[pidInd] < l_MaxProb[retInd-1]) return -1;
	
  //check nsigma cuts
  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)pidInd))>3.0) return -1;
  if(bIsTOFok && l_TOFUsed) if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)pidInd))>3.0) return -1;

  return retInd; //retInd = 1 --> pion, retInd = 2 --> kaon, retInd = 3 --> proton 
}

 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskEffCalc_MC::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
