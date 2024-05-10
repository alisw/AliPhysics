#include "TRandom3.h"
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
#include "AliAnalysisTaskDiffPtFluc_PIDhadrons_v2.h"

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

class AliAnalysisTaskDiffPtFluc_PIDhadrons_v2;
ClassImp(AliAnalysisTaskDiffPtFluc_PIDhadrons_v2)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::AliAnalysisTaskDiffPtFluc_PIDhadrons_v2():
  AliAnalysisTaskSE(),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fInputEvent(0),
  fPIDResponse(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeEvent(0),
  fESDtrackCuts(0),
  fESDtrackCuts_primary(0),
  //fTrigger(AliVEvent::kINT7),
  fMultLow(0),
  fMultHigh(100),
  hNumberOfEvents(0),
  hNumberOfKaonEtaLess0(0),
  hNumberOfPionEtaLess0(0),
  hNumberOfProtonEtaLess0(0),
  fTreeVariableCentrality(0),
  fPtsum_hadrons_less0(0),
  fPtsum_hadrons_greaterEtaMin(0),
  fNsum_hadrons_less0(0),
  fNsum_hadrons_greaterEtaMin(0),
  fNsum_pions_less0(0),
  fNsum_kaons_less0(0),
  fNsum_protons_less0(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  fHistMCEffHadronPlus(0),
  fHistMCEffHadronMinus(0),
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
  f2Dhist_beforeCut_nSigmaTPC_pion(0),
  f2Dhist_beforeCut_nSigmaTPC_kaon(0),
  f2Dhist_beforeCut_nSigmaTPC_proton(0),
  f2Dhist_beforeCut_nSigmaTOF_pion(0),
  f2Dhist_beforeCut_nSigmaTOF_kaon(0),
  f2Dhist_beforeCut_nSigmaTOF_proton(0),
  f2Dhist_beforeCut_nSigmaTPCplusTOF_pion(0),
  f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_beforeCut_nSigmaTPCplusTOF_proton(0),
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion(0),
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon(0),
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton(0),
  f2Dhist_afterCut_nSigmaTPC_pion(0),
  f2Dhist_afterCut_nSigmaTPC_kaon(0),
  f2Dhist_afterCut_nSigmaTPC_proton(0),
  f2Dhist_afterCut_nSigmaTOF_pion(0),
  f2Dhist_afterCut_nSigmaTOF_kaon(0),
  f2Dhist_afterCut_nSigmaTOF_proton(0),
  f2Dhist_afterCut_nSigmaTPCplusTOF_pion(0),
  f2Dhist_afterCut_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_afterCut_nSigmaTPCplusTOF_proton(0),
  f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion(0),
  f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon(0),
  f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton(0),
  f2Dhist_beforeCut_TPCdEdx_all(0),
  f2Dhist_beforeCut_TOFtime_all(0),
  f2Dhist_afterCut_TPCdEdx_all(0),
  f2Dhist_afterCut_TOFtime_all(0),
  f2Dhist_afterCut_TPCdEdx_pion(0),
  f2Dhist_afterCut_TOFtime_pion(0),
  f2Dhist_afterCut_TPCdEdx_kaon(0),
  f2Dhist_afterCut_TOFtime_kaon(0),
  f2Dhist_afterCut_TPCdEdx_proton(0),
  f2Dhist_afterCut_TOFtime_proton(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaCut(0),
  fTPCcrossedrows(0),
  fCentralityEstimator_flag(0),
  fPileupCutVal(0),
  fEtaLeftCut(0),
  fEtaMin(0),
  fEffFlag(0),
  fTreeName(0),
  fEffCorrectionFlag(0),
  fExclusivePIDCut_flag(0),
  fRejectElectron_cut(0),
  fFillTrackQAhists_flag(0),
  fFillPIDhists_flag(0),
  fGlobalTracksAOD(0)
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
   for(int i=0; i<14; i++)
    {
      fPt_no_hadron[i] = 0;
      fPt_no_pion[i] = 0;
      fPt_no_kaon[i] = 0;
      fPt_no_proton[i] = 0;
    }
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::AliAnalysisTaskDiffPtFluc_PIDhadrons_v2(const char *name):
  AliAnalysisTaskSE(name),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fInputEvent(0),
  fPIDResponse(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeEvent(0),
  fESDtrackCuts(0),
  fESDtrackCuts_primary(0),
  //fTrigger(AliVEvent::kINT7),
  fMultLow(0),
  fMultHigh(100),
  hNumberOfEvents(0),
  hNumberOfKaonEtaLess0(0),
  hNumberOfPionEtaLess0(0),
  hNumberOfProtonEtaLess0(0),
  fTreeVariableCentrality(0),
  fPtsum_hadrons_less0(0),
  fPtsum_hadrons_greaterEtaMin(0),
  fNsum_hadrons_less0(0),
  fNsum_hadrons_greaterEtaMin(0),
  fNsum_pions_less0(0),
  fNsum_kaons_less0(0),
  fNsum_protons_less0(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  fHistMCEffHadronPlus(0),
  fHistMCEffHadronMinus(0),
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
  f2Dhist_beforeCut_nSigmaTPC_pion(0),
  f2Dhist_beforeCut_nSigmaTPC_kaon(0),
  f2Dhist_beforeCut_nSigmaTPC_proton(0),
  f2Dhist_beforeCut_nSigmaTOF_pion(0),
  f2Dhist_beforeCut_nSigmaTOF_kaon(0),
  f2Dhist_beforeCut_nSigmaTOF_proton(0),
  f2Dhist_beforeCut_nSigmaTPCplusTOF_pion(0),
  f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_beforeCut_nSigmaTPCplusTOF_proton(0),
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion(0),
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon(0),
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton(0),
  f2Dhist_afterCut_nSigmaTPC_pion(0),
  f2Dhist_afterCut_nSigmaTPC_kaon(0),
  f2Dhist_afterCut_nSigmaTPC_proton(0),
  f2Dhist_afterCut_nSigmaTOF_pion(0),
  f2Dhist_afterCut_nSigmaTOF_kaon(0),
  f2Dhist_afterCut_nSigmaTOF_proton(0),
  f2Dhist_afterCut_nSigmaTPCplusTOF_pion(0),
  f2Dhist_afterCut_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_afterCut_nSigmaTPCplusTOF_proton(0),
  f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion(0),
  f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon(0),
  f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton(0),
  f2Dhist_beforeCut_TPCdEdx_all(0),
  f2Dhist_beforeCut_TOFtime_all(0),
  f2Dhist_afterCut_TPCdEdx_all(0),
  f2Dhist_afterCut_TOFtime_all(0),
  f2Dhist_afterCut_TPCdEdx_pion(0),
  f2Dhist_afterCut_TOFtime_pion(0),
  f2Dhist_afterCut_TPCdEdx_kaon(0),
  f2Dhist_afterCut_TOFtime_kaon(0),
  f2Dhist_afterCut_TPCdEdx_proton(0),
  f2Dhist_afterCut_TOFtime_proton(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaCut(0),
  fTPCcrossedrows(0),
  fCentralityEstimator_flag(0),
  fPileupCutVal(0),
  fEtaLeftCut(0),
  fEtaMin(0),
  fEffFlag(0),
  fTreeName(0),
  fEffCorrectionFlag(0),
  fExclusivePIDCut_flag(0),
  fRejectElectron_cut(0),
  fFillTrackQAhists_flag(0),
  fFillPIDhists_flag(0),
  fGlobalTracksAOD(0)
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
   for(int i=0; i<14; i++)
    {
      fPt_no_hadron[i] = 0;
      fPt_no_pion[i] = 0;
      fPt_no_kaon[i] = 0;
      fPt_no_proton[i] = 0;
    }
  
  fUtils = new AliAnalysisUtils();
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::~AliAnalysisTaskDiffPtFluc_PIDhadrons_v2()  {

  if (fOutputList){
    delete fOutputList;
    fOutputList = 0x0;
  }
  if (fQAList){
    delete fQAList;
    fQAList = 0x0;
  }
  if (fTreeEvent){
    delete fTreeEvent;
    fTreeEvent = 0x0;
  }
  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }
  if(fFBNo==128)
    {
      if(fGlobalTracksAOD) {
	delete fGlobalTracksAOD;
	fGlobalTracksAOD=0x0;
      }
    }
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner(kTRUE);
    fQAList     -> SetOwner(kTRUE);

    OpenFile(1);
    OpenFile(2);
    OpenFile(3);

    
    //QA Plots of Event Selection
    fAODeventCuts.AddQAplotsToList(fQAList,kTRUE);
    fAODeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE, fPileupCutVal);
   
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    // hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    //Number of Kaon finally getting selected with specified cuts event wise
    hNumberOfKaonEtaLess0     = new TH1D ("hNumberOfKaonEtaLess0","",3000,0,3000);
    fOutputList -> Add(hNumberOfKaonEtaLess0);
    //Number of Pion finally getting selected with specified cuts event wise
    hNumberOfPionEtaLess0     = new TH1D ("hNumberOfPionEtaLess0","",3000,0,3000);
    fOutputList -> Add(hNumberOfPionEtaLess0);
    //Number of Proton finally getting selected with specified cuts event wise
    hNumberOfProtonEtaLess0     = new TH1D ("hNumberOfProtonEtaLess0","",3000,0,3000);
    fOutputList -> Add(hNumberOfProtonEtaLess0);

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Histograms for track variables

    //before track cut
    hist_beforeCut_DCAxy = new TH1D("hist_beforeCut_DCAxy","hist_beforeCut_DCAxy", 500, 0, +5);
    hist_beforeCut_DCAz = new TH1D("hist_beforeCut_DCAz","hist_beforeCut_DCAz", 500, 0, +5);
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
    hist_afterCut_DCAxy = new TH1D("hist_afterCut_DCAxy","hist_afterCut_DCAxy", 500, 0, +5);
    hist_afterCut_DCAz = new TH1D("hist_afterCut_DCAz","hist_afterCut_DCAz", 500, 0, +5);
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

    f2Dhist_beforeCut_nSigmaTPC_pion = new TH2D("f2Dhist_beforeCut_nSigmaTPC_pion", "f2Dhist_beforeCut_nSigmaTPC_pion", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPC_kaon = new TH2D("f2Dhist_beforeCut_nSigmaTPC_kaon", "f2Dhist_beforeCut_nSigmaTPC_kaon", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPC_proton = new TH2D("f2Dhist_beforeCut_nSigmaTPC_proton", "f2Dhist_beforeCut_nSigmaTPC_proton", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTOF_pion = new TH2D("f2Dhist_beforeCut_nSigmaTOF_pion", "f2Dhist_beforeCut_nSigmaTOF_pion", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTOF_kaon = new TH2D("f2Dhist_beforeCut_nSigmaTOF_kaon", "f2Dhist_beforeCut_nSigmaTOF_kaon", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTOF_proton = new TH2D("f2Dhist_beforeCut_nSigmaTOF_proton", "f2Dhist_beforeCut_nSigmaTOF_proton", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPCplusTOF_pion = new TH2D("f2Dhist_beforeCut_nSigmaTPCplusTOF_pion", "f2Dhist_beforeCut_nSigmaTPCplusTOF_pion", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon = new TH2D("f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon", "f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon", 500, 0, 5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPCplusTOF_proton = new TH2D("f2Dhist_beforeCut_nSigmaTPCplusTOF_proton", "f2Dhist_beforeCut_nSigmaTPCplusTOF_proton", 500, 0, 5, 2000, -5, +5);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPC_pion);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPC_kaon);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPC_proton);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTOF_pion);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTOF_kaon);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTOF_proton);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPCplusTOF_pion);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPCplusTOF_proton);

    f2Dhist_afterCut_nSigmaTPC_pion = new TH2D("f2Dhist_afterCut_nSigmaTPC_pion", "f2Dhist_afterCut_nSigmaTPC_pion", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPC_kaon = new TH2D("f2Dhist_afterCut_nSigmaTPC_kaon", "f2Dhist_afterCut_nSigmaTPC_kaon", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPC_proton = new TH2D("f2Dhist_afterCut_nSigmaTPC_proton", "f2Dhist_afterCut_nSigmaTPC_proton", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTOF_pion = new TH2D("f2Dhist_afterCut_nSigmaTOF_pion", "f2Dhist_afterCut_nSigmaTOF_pion", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTOF_kaon = new TH2D("f2Dhist_afterCut_nSigmaTOF_kaon", "f2Dhist_afterCut_nSigmaTOF_kaon", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTOF_proton = new TH2D("f2Dhist_afterCut_nSigmaTOF_proton", "f2Dhist_afterCut_nSigmaTOF_proton", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPCplusTOF_pion = new TH2D("f2Dhist_afterCut_nSigmaTPCplusTOF_pion", "f2Dhist_afterCut_nSigmaTPCplusTOF_pion", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPCplusTOF_kaon = new TH2D("f2Dhist_afterCut_nSigmaTPCplusTOF_kaon", "f2Dhist_afterCut_nSigmaTPCplusTOF_kaon", 500, 0, 5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPCplusTOF_proton = new TH2D("f2Dhist_afterCut_nSigmaTPCplusTOF_proton", "f2Dhist_afterCut_nSigmaTPCplusTOF_proton", 500, 0, 5, 2000, -5, +5);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPC_pion);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPC_kaon);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPC_proton);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTOF_pion);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTOF_kaon);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTOF_proton);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPCplusTOF_pion);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPCplusTOF_kaon);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPCplusTOF_proton);

    f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion = new TH2D("f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion", "f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion", 2000, -5, +5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon = new TH2D("f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon", "f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon", 2000, -5, +5, 2000, -5, +5);
    f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton = new TH2D("f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton", "f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton", 2000, -5, +5, 2000, -5, +5);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon);
    fOutputList->Add(f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton);

    f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion = new TH2D("f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion", "f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion", 2000, -5, +5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon = new TH2D("f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon", "f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon", 2000, -5, +5, 2000, -5, +5);
    f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton = new TH2D("f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton", "f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton", 2000, -5, +5, 2000, -5, +5);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon);
    fOutputList->Add(f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton);

    f2Dhist_beforeCut_TPCdEdx_all = new TH2D("f2Dhist_beforeCut_TPCdEdx_all", "f2Dhist_beforeCut_TPCdEdx_all", 500, 0, 5, 10000, 0, 1000);
    f2Dhist_beforeCut_TOFtime_all = new TH2D("f2Dhist_beforeCut_TOFtime_all", "f2Dhist_beforeCut_TOFtime_all", 500, 0, 5, 12000, -6000, 6000);
    
    f2Dhist_afterCut_TPCdEdx_all = new TH2D("f2Dhist_afterCut_TPCdEdx_all", "f2Dhist_afterCut_TPCdEdx_all", 500, 0, 5, 10000, 0, 1000);
    f2Dhist_afterCut_TOFtime_all = new TH2D("f2Dhist_afterCut_TOFtime_all", "f2Dhist_afterCut_TOFtime_all", 500, 0, 5, 12000, -6000, 6000);
    f2Dhist_afterCut_TPCdEdx_pion = new TH2D("f2Dhist_afterCut_TPCdEdx_pion", "f2Dhist_afterCut_TPCdEdx_pion", 500, 0, 5, 10000, 0, 1000);
    f2Dhist_afterCut_TOFtime_pion = new TH2D("f2Dhist_afterCut_TOFtime_pion", "f2Dhist_afterCut_TOFtime_pion", 500, 0, 5, 12000, -6000, 6000);
    f2Dhist_afterCut_TPCdEdx_kaon = new TH2D("f2Dhist_afterCut_TPCdEdx_kaon", "f2Dhist_afterCut_TPCdEdx_kaon", 500, 0, 5, 10000, 0, 1000);
    f2Dhist_afterCut_TOFtime_kaon = new TH2D("f2Dhist_afterCut_TOFtime_kaon", "f2Dhist_afterCut_TOFtime_kaon", 500, 0, 5, 12000, -6000, 6000);
    f2Dhist_afterCut_TPCdEdx_proton = new TH2D("f2Dhist_afterCut_TPCdEdx_proton", "f2Dhist_afterCut_TPCdEdx_proton", 500, 0, 5, 10000, 0, 1000);
    f2Dhist_afterCut_TOFtime_proton = new TH2D("f2Dhist_afterCut_TOFtime_proton", "f2Dhist_afterCut_TOFtime_proton", 500, 0, 5, 12000, -6000, 6000);
    fOutputList->Add(f2Dhist_beforeCut_TPCdEdx_all);
    fOutputList->Add(f2Dhist_beforeCut_TOFtime_all);
    fOutputList->Add(f2Dhist_afterCut_TPCdEdx_all);
    fOutputList->Add(f2Dhist_afterCut_TOFtime_all);
    fOutputList->Add(f2Dhist_afterCut_TPCdEdx_pion);
    fOutputList->Add(f2Dhist_afterCut_TOFtime_pion);
    fOutputList->Add(f2Dhist_afterCut_TPCdEdx_kaon);
    fOutputList->Add(f2Dhist_afterCut_TOFtime_kaon);
    fOutputList->Add(f2Dhist_afterCut_TPCdEdx_proton);
    fOutputList->Add(f2Dhist_afterCut_TOFtime_proton);
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //fTreeEvent: Analysis tree
    
    //TTree object to store variables
    fTreeEvent = new TTree(fTreeName,"Event Tree");
    fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    fTreeEvent->Branch("fPtsum_hadrons_less0",&fPtsum_hadrons_less0,"fPtsum_hadrons_less0/F");
    fTreeEvent->Branch("fPtsum_hadrons_greaterEtaMin",&fPtsum_hadrons_greaterEtaMin,"fPtsum_hadrons_greaterEtaMin/F");
    fTreeEvent->Branch("fNsum_hadrons_less0",&fNsum_hadrons_less0,"fNsum_hadrons_less0/F");
    fTreeEvent->Branch("fNsum_hadrons_greaterEtaMin",&fNsum_hadrons_greaterEtaMin,"fNsum_hadrons_greaterEtaMin/F");
    fTreeEvent->Branch("fNsum_pions_less0",&fNsum_pions_less0,"fNsum_pions_less0/F");
    fTreeEvent->Branch("fNsum_kaons_less0",&fNsum_kaons_less0,"fNsum_kaons_less0/F");
    fTreeEvent->Branch("fNsum_protons_less0",&fNsum_protons_less0,"fNsum_protons_less0/F");
    fTreeEvent->Branch("fPt_no_hadron",&fPt_no_hadron,"fPt_no_hadron[14]/F");
    fTreeEvent->Branch("fPt_no_pion",&fPt_no_pion,"fPt_no_pion[14]/F");
    fTreeEvent->Branch("fPt_no_kaon",&fPt_no_kaon,"fPt_no_kaon[14]/F");
    fTreeEvent->Branch("fPt_no_proton",&fPt_no_proton,"fPt_no_proton[14]/F");

    //----------------------------------------------
    // Look up table for PID information when using TPC-only tracks
    if(fFBNo==128)
      {
	fGlobalTracksAOD = new TExMap();
      }
 
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::UserExec(Option_t *)  {
  
    //Get Input Event
    if ( !GetEvent ()) return;
    //cout<<"*********************** Found AOD event !!! ******************************"<<endl;

    
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
      }

    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();

   
    //Initialize number 0f K+, K-, pi+, pi-, p and p-bar per event
    
 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Int_t ptBinNo = 0;
    Double_t BinCont = 0;
    Double_t EffWgt = 0;

    if(fListTRKCorr)
      {
	//cout<<"## Got efficiency histograms..## "<<endl;
	GetMCEffCorrectionHist();
      }
    else
      cout<<"################ No histograms found ############### "<<endl;



    TH1D *fPt_profile = new TH1D("fPt_profile","fPt_profile", 14, 0.2, 3.0);
    Double_t pT_sum_etaLess0 = 0.0;
    Double_t N_sum_etaLess0 = 0.0;
    Double_t pT_sum_etaGreaterEtamin = 0.0;
    Double_t N_sum_etaGreaterEtamin = 0.0;


    TH1D *fPt_profile_pion = new TH1D("fPt_profile_pion","fPt_profile_pion", 14, 0.2, 3.0);
    TH1D *fPt_profile_kaon = new TH1D("fPt_profile_kaon","fPt_profile_kaon", 14, 0.2, 3.0);
    TH1D *fPt_profile_proton = new TH1D("fPt_profile_proton","fPt_profile_proton", 14, 0.2, 3.0);
    Double_t N_sumPion_etaLess0 = 0.0;
    Double_t N_sumKaon_etaLess0 = 0.0;
    Double_t N_sumProton_etaLess0 = 0.0;

    //Function for efficiency
    TF1 *fEff=new TF1("fEff","[0]*TMath::Exp(-pow([1]/x,[2]))",0.2,3.0);
    fEff->SetParameter(0,0.8);
    fEff->SetParameter(1,0.15);
    fEff->SetParameter(2,1.7);

    Double_t eff, x;
    

    //random no
    TRandom3 ran;

    //Filtering out global tracks for dealing with TPC only tracks
    if(fFBNo==128)
      {
	this->GlobalTracksAOD(fAODevent); 
	if(0 == fGlobalTracksAOD->GetSize())
	  return;
      }


    //Loop on reconstructed tracks
    
    for(Int_t itr=0; itr < fAODevent->GetNumberOfTracks(); itr++)
      {
	
	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;

	Double_t trkPt = aodtrack->Pt();
	Double_t trkPhi = aodtrack->Phi();
	Double_t trkEta = aodtrack->Eta();
	Double_t trkCharge = aodtrack->Charge();
	Double_t trkTPCNCls = aodtrack->GetTPCNcls();
	Double_t trkChi2PerNDF = aodtrack->Chi2perNDF();
	Double_t trkITSchi2 = aodtrack->GetITSchi2();
	Int_t trkITSNcls = aodtrack->GetITSNcls();
	Double_t trkITSchi2perNcls = trkITSchi2/trkITSNcls;
	Double_t trkTPCchi2perNcls = aodtrack->GetTPCchi2perCluster();
	Double_t trkTPCcrossedrows = aodtrack->GetTPCCrossedRows();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Double_t trkDCAxy = aodtrack->DCA();
	// Double_t trkDCAz = aodtrack->ZAtDCA();
	
	float vertexX = -999.;
	float vertexY = -999.;
	float vertexZ = -999.;

	if(vertex) {
	  Double32_t fCov[6];
	  vertex->GetCovarianceMatrix(fCov);
	  if(vertex->GetNContributors() > 0) {
	    if(fCov[5] != 0) {
	      vertexX = vertex->GetX();
	      vertexY = vertex->GetY();
	      vertexZ = vertex->GetZ();
	    }
	  }
	}
	Double_t pos[3];
	aodtrack->GetXYZ(pos);
	Double_t DCAX, DCAY, DCAZ, trkDCAxy, trkDCAz;
	DCAX = pos[0] - vertexX;
	DCAY = pos[1] - vertexY;
	DCAZ = pos[2] - vertexZ;
	trkDCAxy = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));
	trkDCAz = TMath::Abs(DCAZ);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	//Filling track QA before trackselection cuts
	if(fFillTrackQAhists_flag==1)
	  {
	    FilltrackQAplots_beforeCut(trkDCAxy, trkDCAz, trkEta, trkITSchi2perNcls, trkTPCchi2perNcls, trkTPCcrossedrows);
	  }
	
	
	//TRACK SELECTION CUTS

	//Track selectionL FilterBit :   default-96, systematics-768
	if(!aodtrack->TestFilterBit(fFBNo))  continue;
	//if(!aodtrack->TestFilterBit(768))  continue;
	
	//cuts on TPCchi2perClstr and ITSchi2perClstr and TPC nCrossedRows
	if (trkTPCcrossedrows < fTPCcrossedrows)
	  continue;
	if (trkTPCchi2perNcls > fChi2TPC)
	  continue;
	if (trkITSchi2perNcls > fChi2ITS)
	  continue;
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Filling PID QA plots before applying PID cuts
	
	if(fFillPIDhists_flag==1)
	  {
	    FillPIDQAplots_beforeCut(track);
	  }
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Kinematic cuts on pT and Eta
	if (TMath::Abs(trkEta) > 0.8) continue;
	if (trkPt < 0.2) continue;
	if (trkPt > 3.0) continue;


	//Filling track QA before trackselection cuts
	if(fFillTrackQAhists_flag==1)
	  {
	    FilltrackQAplots_afterCut(trkDCAxy, trkDCAz, trkEta, trkITSchi2perNcls, trkTPCchi2perNcls, trkTPCcrossedrows);
	  }
	
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++
	int flag=1;
	//Particle loss imposing

	//using algebraic function
	if(fEffFlag==1)
	  {
	    x=ran.Uniform(0,1);
	    eff=fEff->Eval(trkPt);
	    //cout<<x<<"\t"<<eff<<endl;
	    if(x > eff)
	      flag = 0;
	  }
	//using actual efficiencies from files
	if(fEffFlag==2)
	  {
	    if(trkCharge < 0)
	      {
		ptBinNo = fHistMCEffHadronMinus->FindBin(trkPt);
		BinCont = fHistMCEffHadronMinus->GetBinContent(ptBinNo);
	      }
	    if(trkCharge > 0)
	      {
		ptBinNo = fHistMCEffHadronPlus->FindBin(trkPt);
		BinCont = fHistMCEffHadronPlus->GetBinContent(ptBinNo);
	      }
	    x=ran.Uniform(0,1);
	    if(x > BinCont)
	      flag = 0;
	  }
	
	if(flag==0)
	  continue;
	//+++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++

	
	if(fEffCorrectionFlag == 0)
	  {
	    if(TMath::Abs(trkCharge) > 0)
	      {
		if(trkEta < fEtaLeftCut)
		  {
		    fPt_profile->Fill(trkPt);
		    pT_sum_etaLess0 += trkPt;
		    N_sum_etaLess0 += 1.0;
		  }
		if(trkEta > fEtaMin) //fEtaMin is right boundary of EtaGap
		  {
		    pT_sum_etaGreaterEtamin += trkPt;
		    N_sum_etaGreaterEtamin += 1.0;
		  }
	      }
	  }

	if(fEffCorrectionFlag == 1)
	  {
	    if(TMath::Abs(trkCharge) > 0)
	      {
		if(trkCharge < 0)
		  {
		    ptBinNo = fHistMCEffHadronMinus->FindBin(trkPt);
		    BinCont = fHistMCEffHadronMinus->GetBinContent(ptBinNo);
		    //cout<<"Eff: "<<BinCont<<endl;
		  }
		if(trkCharge > 0)
		  {
		    ptBinNo = fHistMCEffHadronPlus->FindBin(trkPt);
		    BinCont = fHistMCEffHadronPlus->GetBinContent(ptBinNo);
		  }
		
		if(trkEta < fEtaLeftCut)
		  {
		    if(BinCont != 0)
		      {
			fPt_profile->Fill(trkPt,1.0/BinCont);
			pT_sum_etaLess0 += trkPt/BinCont;
			N_sum_etaLess0 += 1.0/BinCont;
		      }
		  }
	    
		if(trkEta > fEtaMin)
		  {
		    if(BinCont != 0)
		      {
			pT_sum_etaGreaterEtamin += trkPt/BinCont;
			N_sum_etaGreaterEtamin += 1.0/BinCont;
		      }
		  }
	      }
	  }

	if(fFBNo==128)
	  {
	    Int_t id=track->GetID();
	    Int_t retrievekey=-id-1;
	    Int_t iGlobalTrackIndex=fGlobalTracksAOD->GetValue(retrievekey);
	    track = (AliVTrack*)fAODevent->GetTrack(iGlobalTrackIndex);
	    if(!track)      continue;
	    aodtrack  = dynamic_cast<AliAODTrack*>(track);
	    if(!aodtrack)      continue;
	  }

	
	//PID selection
	Bool_t IsKaon = KaonSelector(track, fPIDnSigmaCut);
	Bool_t IsPion = PionSelector(track, fPIDnSigmaCut);
	Bool_t IsProton = ProtonSelector(track, fPIDnSigmaCut);

	if (!IsKaon && !IsPion && !IsProton) continue;

	if(fExclusivePIDCut_flag==1)
	  {
	    //Check if a particle is selected as more than one type
	    Int_t flag = 0;
	    if(IsPion) flag+=1;
	    if(IsKaon) flag+=1;
	    if(IsProton) flag+=1;
	    if(flag>1)
	      {
		cout<<"Particle identified as more than on PID: flag= "<<flag<<endl;
		continue;
	      }
	  }

	if(fRejectElectron_cut > 0)
	  {
	    Bool_t IsElectron = ElectronRejectionCut(track, fRejectElectron_cut);
	    if(IsElectron) continue;
	  }

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Filling PID QA plots after applying PID cuts
	
	if(fFillPIDhists_flag==1)
	  {
	    FillPIDQAplots_afterCut(track, IsPion, IsKaon, IsProton);
	  }
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if(fEffCorrectionFlag == 0)
	  {
	    if(TMath::Abs(trkCharge) > 0)
	      {
		if(trkEta < fEtaLeftCut)
		  {
		    if(IsPion)
		      {
			fPt_profile_pion->Fill(trkPt);
			N_sumPion_etaLess0 += 1.0;
		      }
		    if(IsKaon)
		      {
			fPt_profile_kaon->Fill(trkPt);
			N_sumKaon_etaLess0 += 1.0;
		      }
		    if(IsProton && trkPt > 0.4)
		      {
			fPt_profile_proton->Fill(trkPt);
			N_sumProton_etaLess0 += 1.0;
		      }
		  }
	      }
	  }

	if(fEffCorrectionFlag == 1)
	  {
	    if(TMath::Abs(trkCharge) > 0)
	      {
		if(trkEta < fEtaLeftCut)
		  {
		    if(IsPion)
		      {
			if(trkCharge < 0)
			  {
			    ptBinNo = fHistMCEffPionMinus->FindBin(trkPt);
			    BinCont = fHistMCEffPionMinus->GetBinContent(ptBinNo);
			  }
			if(trkCharge > 0)
			  {
			    ptBinNo = fHistMCEffPionPlus->FindBin(trkPt);
			    BinCont = fHistMCEffPionPlus->GetBinContent(ptBinNo);
			  }
			if(BinCont != 0)
			  {
			    fPt_profile_pion->Fill(trkPt,1.0/BinCont);
			    N_sumPion_etaLess0 += 1.0/BinCont;
			  }
		      }
		    
		    if(IsKaon)
		      {
			if(trkCharge < 0)
			  {
			    ptBinNo = fHistMCEffKaonMinus->FindBin(trkPt);
			    BinCont = fHistMCEffKaonMinus->GetBinContent(ptBinNo);
			  }
			if(trkCharge > 0)
			  {
			    ptBinNo = fHistMCEffKaonPlus->FindBin(trkPt);
			    BinCont = fHistMCEffKaonPlus->GetBinContent(ptBinNo);
			  }
			if(BinCont != 0)
			  {
			    fPt_profile_kaon->Fill(trkPt,1.0/BinCont);
			    N_sumKaon_etaLess0 += 1.0/BinCont;
			  }
		      }
		    
		    if(IsProton && trkPt > 0.4)
		      {
			if(trkCharge < 0)
			  {
			    ptBinNo = fHistMCEffProtonMinus->FindBin(trkPt);
			    BinCont = fHistMCEffProtonMinus->GetBinContent(ptBinNo);
			  }
			if(trkCharge > 0)
			  {
			    ptBinNo = fHistMCEffProtonPlus->FindBin(trkPt);
			    BinCont = fHistMCEffProtonPlus->GetBinContent(ptBinNo);
			  }
			if(BinCont != 0)
			  {
			    fPt_profile_proton->Fill(trkPt,1.0/BinCont);
			    N_sumProton_etaLess0 += 1.0/BinCont;
			  }
		      }
		  }
	      }
	  }
	
      }      
    //end reconstructed track loop
    
    
    
    //Kaon
    hNumberOfKaonEtaLess0->Fill(N_sumKaon_etaLess0);
    //Pion
    hNumberOfPionEtaLess0->Fill(N_sumPion_etaLess0);
    //Proton
    hNumberOfProtonEtaLess0->Fill(N_sumProton_etaLess0);

    //Tree Variables++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    fTreeVariableCentrality=lV0M;
    fPtsum_hadrons_less0=pT_sum_etaLess0;
    fPtsum_hadrons_greaterEtaMin=pT_sum_etaGreaterEtamin;
    fNsum_hadrons_less0=N_sum_etaLess0;
    fNsum_hadrons_greaterEtaMin=N_sum_etaGreaterEtamin;
    fNsum_pions_less0=N_sumPion_etaLess0;
    fNsum_kaons_less0=N_sumKaon_etaLess0;
    fNsum_protons_less0=N_sumProton_etaLess0;
    
    for(int i=0; i<14; i++)
      {
	fPt_no_hadron[i]=fPt_profile->GetBinContent(i+1);

	fPt_no_pion[i]=fPt_profile_pion->GetBinContent(i+1);
	fPt_no_kaon[i]=fPt_profile_kaon->GetBinContent(i+1);
	fPt_no_proton[i]=fPt_profile_proton->GetBinContent(i+1);
      }
    
    fTreeEvent->Fill();

    fPt_profile->Delete();
    fPt_profile_pion->Delete();
    fPt_profile_kaon->Delete();
    fPt_profile_proton->Delete();
  
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::GetEvent ()  //event cuts copied from my code written earlier 

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
	  PostData(3, fTreeEvent);
	  return kFALSE;
	}
    }

  fInputEvent = dynamic_cast <AliVEvent*>(InputEvent());
  if(!fInputEvent)
    {
      AliWarning("ERROR: fInputEvent (AliVEvent) not available \n");
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  
  hNumberOfEvents -> Fill(0.5);

    
  //Standard Event Cuts
  if (!fAODeventCuts.AcceptEvent(fInputEvent)) {
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
    return kFALSE;
  }
  hNumberOfEvents -> Fill(1.5);
  
  
  //Reject Events with Incomplete DAQ
  if (fAODevent->IsIncompleteDAQ())
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(2.5);
        
  //V0 Timing Decision
  AliVVZERO *vzeroData = fInputEvent->GetVZEROData();
  if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision()))
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
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
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(4.5);

    
  //Primary Vertex Tracks
  AliAODVertex *vertex_tracks = (AliAODVertex*) fAODevent->GetPrimaryVertexTracks();
  if (!vertex_tracks)
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(5.5);
  
        
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(6.5);

    
        
    //Primary Vertex SPD
    AliAODVertex *vertex_SPD = (AliAODVertex*) fAODevent->GetPrimaryVertexSPD();
    if (!vertex_SPD)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(7.5);
        
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(8.5);

    
    //SPD Pile-up in Mult Bins
    if (fAODevent->IsPileupFromSPDInMultBins())
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
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
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(10.5);
    
        
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(11.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -fVertexZMax || vertex_tracks->GetZ() > +fVertexZMax)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(12.5);
               
    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      } 
    hNumberOfEvents -> Fill(13.5);
                
    //Selection of Multiplicity Range
    Double_t mult_percentile = -999.0;
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
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(14.5);
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    return kTRUE;
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::PassedTrackQualityCuts (AliAODTrack *track)  {
    
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::KaonSelector(AliVTrack *track, Double_t nSigmaCut)  {
 
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
      if(TMath::Abs(fTPCnSigmaKaon) < nSigmaCut)
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
      if (fTPCplusTOFnSigmaKaon < nSigmaCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
  
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::PionSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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
      
      if(TMath::Abs(fTPCnSigmaPion) < nSigmaCut)
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
      
      if (fTPCplusTOFnSigmaPion < nSigmaCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::ProtonSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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

      if(TMath::Abs(fTPCnSigmaProton) < nSigmaCut)
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
      
      if (fTPCplusTOFnSigmaProton < nSigmaCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::ElectronRejectionCut(AliVTrack *track, Int_t fCut)
{
  //TPC nsigma
  Double_t fTPCnSigma_Pion = 0.0;
  Double_t fTPCnSigma_Proton = 0.0;
  Double_t fTPCnSigma_Kaon = 0.0;
  Double_t fTPCnSigma_Electron = 0.0;
  fTPCnSigma_Pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigma_Proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigma_Kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  fTPCnSigma_Electron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

  /*
  //TOF nsigma
  Double_t fTOFnSigma_Pion = 0.0;
  Double_t fTOFnSigma_Proton = 0.0;
  Double_t fTOFnSigma_Kaon = 0.0;
  Double_t fTOFnSigma_Electron = 0.0;	
  fTOFnSigma_Pion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigma_Proton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigma_Kaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
  fTOFnSigma_Electron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
  */

  if(fCut==1) //Cut 1
    {
      if(TMath::Abs(fTPCnSigma_Electron) < 3.0 && TMath::Abs(fTPCnSigma_Pion) > 3.0 && TMath::Abs(fTPCnSigma_Kaon) > 3.0 && TMath::Abs(fTPCnSigma_Proton) > 3.0)
	return kTRUE;
    }
  else if(fCut==2) //Cut 2
    {
      if(TMath::Abs(fTPCnSigma_Electron) < 1.0)
	return kTRUE;
    }
  else
    return kFALSE;

}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kTRUE);
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  return passedTrackPileupCut;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::GetMCEffCorrectionHist()
{
  if(fListTRKCorr)
    {
      fHistMCEffKaonPlus = (TH1D*) fListTRKCorr->FindObject("histKaonPlusEff");
      fHistMCEffKaonMinus = (TH1D*) fListTRKCorr->FindObject("histKaonMinusEff");

      fHistMCEffPionPlus = (TH1D*) fListTRKCorr->FindObject("histPionPlusEff");
      fHistMCEffPionMinus = (TH1D*) fListTRKCorr->FindObject("histPionMinusEff");

      fHistMCEffProtonPlus = (TH1D*) fListTRKCorr->FindObject("histProtonPlusEff");
      fHistMCEffProtonMinus = (TH1D*) fListTRKCorr->FindObject("histProtonMinusEff");

      fHistMCEffHadronPlus = (TH1D*) fListTRKCorr->FindObject("histHadronPlusEff");
      fHistMCEffHadronMinus = (TH1D*) fListTRKCorr->FindObject("histHadronMinusEff");

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
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::FilltrackQAplots_beforeCut(Double_t fDcaXY, Double_t fDcaZ, Double_t fEta, Double_t fITSchi2perNcls, Double_t fTPCchi2perNcls, Double_t fTPCcrossedrows)
{
  hist_beforeCut_DCAxy->Fill(fDcaXY);
  hist_beforeCut_DCAz->Fill(fDcaZ);
  hist_beforeCut_eta->Fill(fEta);
  hist_beforeCut_chi2perTPCclstr->Fill(fTPCchi2perNcls);
  hist_beforeCut_chi2perITSclstr->Fill(fITSchi2perNcls);
  hist_beforeCut_TPCncrossedrows->Fill(fTPCcrossedrows);
}

//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::FilltrackQAplots_afterCut(Double_t fDcaXY, Double_t fDcaZ, Double_t fEta, Double_t fITSchi2perNcls, Double_t fTPCchi2perNcls, Double_t fTPCcrossedrows)
{
  hist_afterCut_DCAxy->Fill(fDcaXY);
  hist_afterCut_DCAz->Fill(fDcaZ);
  hist_afterCut_eta->Fill(fEta);
  hist_afterCut_chi2perTPCclstr->Fill(fTPCchi2perNcls);
  hist_afterCut_chi2perITSclstr->Fill(fITSchi2perNcls);
  hist_afterCut_TPCncrossedrows->Fill(fTPCcrossedrows);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::FillPIDQAplots_beforeCut(AliVTrack *track)
{
  Double_t p[3];
  track->PxPyPz(p);
  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fdEdx = track->GetTPCsignal();
  Double_t fTofSig = track->GetTOFsignal();
  Double_t fPidTime[5];
  track->GetIntegratedTimes(fPidTime);

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

  f2Dhist_beforeCut_TPCdEdx_all->Fill(Track_pt, fdEdx);
  f2Dhist_beforeCut_TOFtime_all->Fill(Track_pt, fTofSig-fPidTime[2]);

  //TPC
  f2Dhist_beforeCut_nSigmaTPC_pion->Fill(Track_pt, fTPCnSigma_Pion);
  f2Dhist_beforeCut_nSigmaTPC_kaon->Fill(Track_pt, fTPCnSigma_Kaon);
  f2Dhist_beforeCut_nSigmaTPC_proton->Fill(Track_pt, fTPCnSigma_Proton);
  //TOF
  f2Dhist_beforeCut_nSigmaTOF_pion->Fill(Track_pt, fTOFnSigma_Pion);
  f2Dhist_beforeCut_nSigmaTOF_kaon->Fill(Track_pt, fTOFnSigma_Kaon);
  f2Dhist_beforeCut_nSigmaTOF_proton->Fill(Track_pt, fTOFnSigma_Proton);
  //TPC+TOF
  f2Dhist_beforeCut_nSigmaTPCplusTOF_pion->Fill(Track_pt, fTPCplusTOFnSigma_Pion);
  f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon->Fill(Track_pt, fTPCplusTOFnSigma_Kaon);
  f2Dhist_beforeCut_nSigmaTPCplusTOF_proton->Fill(Track_pt, fTPCplusTOFnSigma_Proton);
  //TPC vs. TOF
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion->Fill(fTPCnSigma_Pion, fTOFnSigma_Pion);
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon->Fill(fTPCnSigma_Kaon, fTOFnSigma_Kaon);
  f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton->Fill(fTPCnSigma_Proton, fTOFnSigma_Proton);

}

//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::FillPIDQAplots_afterCut(AliVTrack *track, Bool_t Pion_flag, Bool_t Kaon_flag, Bool_t Proton_flag)
{
  Double_t p[3];
  track->PxPyPz(p);
  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fdEdx = track->GetTPCsignal();
  Double_t fTofSig = track->GetTOFsignal();
  Double_t fPidTime[5];
  track->GetIntegratedTimes(fPidTime);

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

  f2Dhist_afterCut_TPCdEdx_all->Fill(Track_pt, fdEdx);
  f2Dhist_afterCut_TOFtime_all->Fill(Track_pt, fTofSig-fPidTime[2]);
	
  if(Pion_flag)
    {
      f2Dhist_afterCut_TPCdEdx_pion->Fill(Track_pt, fdEdx);
      f2Dhist_afterCut_TOFtime_pion->Fill(Track_pt, fTofSig-fPidTime[2]);
      f2Dhist_afterCut_nSigmaTPC_pion->Fill(Track_pt, fTPCnSigma_Pion);
      f2Dhist_afterCut_nSigmaTOF_pion->Fill(Track_pt, fTOFnSigma_Pion);
      f2Dhist_afterCut_nSigmaTPCplusTOF_pion->Fill(Track_pt, fTPCplusTOFnSigma_Pion);
      f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion->Fill(fTPCnSigma_Pion, fTOFnSigma_Pion);
    }

  if(Kaon_flag)
    {
      f2Dhist_afterCut_TPCdEdx_kaon->Fill(Track_pt, fdEdx);
      f2Dhist_afterCut_TOFtime_kaon->Fill(Track_pt, fTofSig-fPidTime[3]);
      f2Dhist_afterCut_nSigmaTPC_kaon->Fill(Track_pt, fTPCnSigma_Kaon);
      f2Dhist_afterCut_nSigmaTOF_kaon->Fill(Track_pt, fTOFnSigma_Kaon);
      f2Dhist_afterCut_nSigmaTPCplusTOF_kaon->Fill(Track_pt, fTPCplusTOFnSigma_Kaon);
      f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon->Fill(fTPCnSigma_Kaon, fTOFnSigma_Kaon);
    }

  if(Proton_flag)
    {
      f2Dhist_afterCut_TPCdEdx_proton->Fill(Track_pt, fdEdx);
      f2Dhist_afterCut_TOFtime_proton->Fill(Track_pt, fTofSig-fPidTime[4]);
      f2Dhist_afterCut_nSigmaTPC_proton->Fill(Track_pt, fTPCnSigma_Proton);
      f2Dhist_afterCut_nSigmaTOF_proton->Fill(Track_pt, fTOFnSigma_Proton);
      f2Dhist_afterCut_nSigmaTPCplusTOF_proton->Fill(Track_pt, fTPCplusTOFnSigma_Proton);
      f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton->Fill(fTPCnSigma_Proton, fTOFnSigma_Proton);
    }
	
}
//______________________________________________________________________________________________________________________________________ **From Ante**
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::GlobalTracksAOD(AliAODEvent *aAOD)
{
 // Filter out unique global tracks in AOD and store them in fGlobalTracksAOD.

 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex.
 //           However, only the "normal" global tracks come with positive ID, additionally they can be discriminated simply via: aodTrack->IsGlobalConstrained()
 //           Global constrained tracks have the same negative ID as the TPC-only tracks, both associated with the same "normal global" tracks. E.g. we can have
 //            iTrack: atrack->GetID(): atrack->Pt() atrack->Eta() atrack->Phi()
 //                 1:               0:     2.073798     -0.503640      2.935432
 //                19:              -1:     2.075537     -0.495988      2.935377 => this is TPC-only
 //                35:              -1:     2.073740     -0.493576      2.935515 => this is IsGlobalConstrained()
 //           In fact, this is important, otherwise there is double or even triple counting in some cases.
 // Remark 2: There are tracks for which: 0 == aodTrack->GetFilterMap()
 //           a) Basically all of them pass: atrack->GetType() == AliAODTrack::kFromDecayVtx , but few exceptions also pass atrack->GetType() == AliAODTrack::kPrimary
 //           b) All of them apparently have positive ID, i.e. these are global tracks
 //           c) Clearly, we cannot use TestFilterBit() on them
 //           d) None of them apparently satisfies: atrack->IsGlobalConstrained()
 // Remark 3: There is a performance penalty when fGlobalTracksAOD[1] and fGlobalTracksAOD[2] needed for mixed events are calculated.
 //           Yes, I can get them directly from fGlobalTracksAOD[0], without calling this method for them again. TBI today

 // a) Insanity checks;
 // b) Determine the map.

 // a) Insanity checks:
 if(0 != fGlobalTracksAOD->GetSize()){fGlobalTracksAOD->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:

 //if(fUseFisherYates){cout<<__LINE__<<endl;exit(1);} // TBI 20210810 check and validate if also here Fisher-Yates needs to be applied

 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAOD->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

}
 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons_v2::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
