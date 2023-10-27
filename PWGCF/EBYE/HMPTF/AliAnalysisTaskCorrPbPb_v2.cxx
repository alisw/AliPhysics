
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
#include "AliAnalysisTaskCorrPbPb_v2.h"

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

class AliAnalysisTaskCorrPbPb_v2;
ClassImp(AliAnalysisTaskCorrPbPb_v2)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPbPb_v2::AliAnalysisTaskCorrPbPb_v2():
  AliAnalysisTaskSE(),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fInputEvent(0),
  fPIDResponse(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeList(0),
  //fTreeEvent(0),
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
  // fTreeVariableCentrality(0),
  // fNoKaonPlus_ptmax2(0),
  // fNoKaonMinus_ptmax2(0),
  // fNoKaonPlus_ptmax3(0),
  // fNoKaonMinus_ptmax3(0),
  // fNoPionPlus_ptmax2(0),
  // fNoPionMinus_ptmax2(0),
  // fNoPionPlus_ptmax3(0),
  // fNoPionMinus_ptmax3(0),
  // fNoProtonPlus_ptmax2(0),
  // fNoProtonMinus_ptmax2(0),
  // fNoProtonPlus_ptmax3(0),
  // fNoProtonMinus_ptmax3(0),
  // fCorrectedNoKaonPlus_ptmax2(0),
  // fCorrectedNoKaonMinus_ptmax2(0),
  // fCorrectedNoKaonPlus_ptmax3(0),
  // fCorrectedNoKaonMinus_ptmax3(0),
  // fCorrectedNoPionPlus_ptmax2(0),
  // fCorrectedNoPionMinus_ptmax2(0),
  // fCorrectedNoPionPlus_ptmax3(0),
  // fCorrectedNoPionMinus_ptmax3(0),
  // fCorrectedNoProtonPlus_ptmax2(0),
  // fCorrectedNoProtonMinus_ptmax2(0),
  // fCorrectedNoProtonPlus_ptmax3(0),
  // fCorrectedNoProtonMinus_ptmax3(0),
  // fEffSqrFactrPionMinus_ptmax2(0),
  // fEffSqrFactrPionPlus_ptmax2(0),
  // fEffSqrFactrProtonMinus_ptmax2(0),
  // fEffSqrFactrProtonPlus_ptmax2(0),
  // fEffSqrFactrKaonMinus_ptmax2(0),
  // fEffSqrFactrKaonPlus_ptmax2(0),
  // fEffSqrFactrPionMinus_ptmax3(0),
  // fEffSqrFactrPionPlus_ptmax3(0),
  // fEffSqrFactrProtonMinus_ptmax3(0),
  // fEffSqrFactrProtonPlus_ptmax3(0),
  // fEffSqrFactrKaonMinus_ptmax3(0),
  // fEffSqrFactrKaonPlus_ptmax3(0),
  Profile_ptmax2_RecNetKaon(0),
  Profile_ptmax2_RecNetProton(0),
  Profile_ptmax2_RecNetCharge(0),
  Profile_ptmax2_RecNetChargeNetKaon(0),
  Profile_ptmax2_RecNetChargeNetProton(0),
  Profile_ptmax2_RecNetKaonNetProton(0),
  Profile_ptmax2_Term1C2RecNetKaon(0),
  Profile_ptmax2_Term2C2RecNetKaon(0),
  Profile_ptmax2_Term1C2RecNetCharge(0),
  Profile_ptmax2_Term2C2RecNetCharge(0),
  Profile_ptmax2_Term1C2RecNetProton(0),
  Profile_ptmax2_Term2C2RecNetProton(0),
  Profile_ptSTAR_RecNetKaon(0),
  Profile_ptSTAR_RecNetProton(0),
  Profile_ptSTAR_RecNetCharge(0),
  Profile_ptSTAR_RecNetChargeNetKaon(0),
  Profile_ptSTAR_RecNetChargeNetProton(0),
  Profile_ptSTAR_RecNetKaonNetProton(0),
  Profile_ptSTAR_Term1C2RecNetKaon(0),
  Profile_ptSTAR_Term2C2RecNetKaon(0),
  Profile_ptSTAR_Term1C2RecNetCharge(0),
  Profile_ptSTAR_Term2C2RecNetCharge(0),
  Profile_ptSTAR_Term1C2RecNetProton(0),
  Profile_ptSTAR_Term2C2RecNetProton(0),
  Profile_ptmax2_CorrectedNetKaon(0),
  Profile_ptmax2_CorrectedNetProton(0),
  Profile_ptmax2_CorrectedNetCharge(0),
  Profile_ptmax2_CorrectedNetChargeNetKaon_term1(0),
  Profile_ptmax2_CorrectedNetChargeNetKaon_term2(0),
  Profile_ptmax2_CorrectedNetChargeNetKaon_term3(0),
  Profile_ptmax2_CorrectedNetChargeNetProton_term1(0),
  Profile_ptmax2_CorrectedNetChargeNetProton_term2(0),
  Profile_ptmax2_CorrectedNetChargeNetProton_term3(0),
  Profile_ptmax2_CorrectedNetKaonNetProton(0),
  Profile_ptmax2_Term1C2CorrectedNetKaon(0),
  Profile_ptmax2_Term2C2CorrectedNetKaon(0),
  Profile_ptmax2_Term3C2CorrectedNetKaon(0),
  Profile_ptmax2_Term4C2CorrectedNetKaon(0),
  Profile_ptmax2_Term1C2CorrectedNetCharge(0),
  Profile_ptmax2_Term2C2CorrectedNetCharge(0),
  Profile_ptmax2_Term3C2CorrectedNetCharge(0),
  Profile_ptmax2_Term4C2CorrectedNetCharge(0),
  Profile_ptmax2_Term1C2CorrectedNetProton(0),
  Profile_ptmax2_Term2C2CorrectedNetProton(0),
  Profile_ptmax2_Term3C2CorrectedNetProton(0),
  Profile_ptmax2_Term4C2CorrectedNetProton(0),
  Profile_ptSTAR_CorrectedNetKaon(0),
  Profile_ptSTAR_CorrectedNetProton(0),
  Profile_ptSTAR_CorrectedNetCharge(0),
  Profile_ptSTAR_CorrectedNetChargeNetKaon_term1(0),
  Profile_ptSTAR_CorrectedNetChargeNetKaon_term2(0),
  Profile_ptSTAR_CorrectedNetChargeNetKaon_term3(0),
  Profile_ptSTAR_CorrectedNetChargeNetProton_term1(0),
  Profile_ptSTAR_CorrectedNetChargeNetProton_term2(0),
  Profile_ptSTAR_CorrectedNetChargeNetProton_term3(0),
  Profile_ptSTAR_CorrectedNetKaonNetProton(0),
  Profile_ptSTAR_Term1C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term2C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term3C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term4C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term1C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term2C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term3C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term4C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term1C2CorrectedNetProton(0),
  Profile_ptSTAR_Term2C2CorrectedNetProton(0),
  Profile_ptSTAR_Term3C2CorrectedNetProton(0),
  Profile_ptSTAR_Term4C2CorrectedNetProton(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaPionCut(0),
  fPIDnSigmaKaonCut(0),
  fPIDnSigmaProtonCut(0),
  fTPCcrossedrows(0),
  fEtaMax(0),
  fPileupCutVal(0),
  fCentralityEstimator_flag(0),
  f2Dhist_nSigmaTPC_pion(0),
  f2Dhist_nSigmaTPC_kaon(0),
  f2Dhist_nSigmaTPC_proton(0),
  f2Dhist_nSigmaTOF_pion(0),
  f2Dhist_nSigmaTOF_kaon(0),
  f2Dhist_nSigmaTOF_proton(0),
  f2Dhist_nSigmaTPCplusTOF_pion(0),
  f2Dhist_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_nSigmaTPCplusTOF_proton(0),
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
  hist_afterCut_TPCncrossedrows(0)
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
AliAnalysisTaskCorrPbPb_v2::AliAnalysisTaskCorrPbPb_v2(const char *name):
  AliAnalysisTaskSE(name),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fInputEvent(0),
  fPIDResponse(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeList(0),
  //fTreeEvent(0),
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
  // fTreeVariableCentrality(0),
  // fNoKaonPlus_ptmax2(0),
  // fNoKaonMinus_ptmax2(0),
  // fNoKaonPlus_ptmax3(0),
  // fNoKaonMinus_ptmax3(0),
  // fNoPionPlus_ptmax2(0),
  // fNoPionMinus_ptmax2(0),
  // fNoPionPlus_ptmax3(0),
  // fNoPionMinus_ptmax3(0),
  // fNoProtonPlus_ptmax2(0),
  // fNoProtonMinus_ptmax2(0),
  // fNoProtonPlus_ptmax3(0),
  // fNoProtonMinus_ptmax3(0),
  // fCorrectedNoKaonPlus_ptmax2(0),
  // fCorrectedNoKaonMinus_ptmax2(0),
  // fCorrectedNoKaonPlus_ptmax3(0),
  // fCorrectedNoKaonMinus_ptmax3(0),
  // fCorrectedNoPionPlus_ptmax2(0),
  // fCorrectedNoPionMinus_ptmax2(0),
  // fCorrectedNoPionPlus_ptmax3(0),
  // fCorrectedNoPionMinus_ptmax3(0),
  // fCorrectedNoProtonPlus_ptmax2(0),
  // fCorrectedNoProtonMinus_ptmax2(0),
  // fCorrectedNoProtonPlus_ptmax3(0),
  // fCorrectedNoProtonMinus_ptmax3(0),
  // fEffSqrFactrPionMinus_ptmax2(0),
  // fEffSqrFactrPionPlus_ptmax2(0),
  // fEffSqrFactrProtonMinus_ptmax2(0),
  // fEffSqrFactrProtonPlus_ptmax2(0),
  // fEffSqrFactrKaonMinus_ptmax2(0),
  // fEffSqrFactrKaonPlus_ptmax2(0),
  // fEffSqrFactrPionMinus_ptmax3(0),
  // fEffSqrFactrPionPlus_ptmax3(0),
  // fEffSqrFactrProtonMinus_ptmax3(0),
  // fEffSqrFactrProtonPlus_ptmax3(0),
  // fEffSqrFactrKaonMinus_ptmax3(0),
  // fEffSqrFactrKaonPlus_ptmax3(0),
  Profile_ptmax2_RecNetKaon(0),
  Profile_ptmax2_RecNetProton(0),
  Profile_ptmax2_RecNetCharge(0),
  Profile_ptmax2_RecNetChargeNetKaon(0),
  Profile_ptmax2_RecNetChargeNetProton(0),
  Profile_ptmax2_RecNetKaonNetProton(0),
  Profile_ptmax2_Term1C2RecNetKaon(0),
  Profile_ptmax2_Term2C2RecNetKaon(0),
  Profile_ptmax2_Term1C2RecNetCharge(0),
  Profile_ptmax2_Term2C2RecNetCharge(0),
  Profile_ptmax2_Term1C2RecNetProton(0),
  Profile_ptmax2_Term2C2RecNetProton(0),
  Profile_ptSTAR_RecNetKaon(0),
  Profile_ptSTAR_RecNetProton(0),
  Profile_ptSTAR_RecNetCharge(0),
  Profile_ptSTAR_RecNetChargeNetKaon(0),
  Profile_ptSTAR_RecNetChargeNetProton(0),
  Profile_ptSTAR_RecNetKaonNetProton(0),
  Profile_ptSTAR_Term1C2RecNetKaon(0),
  Profile_ptSTAR_Term2C2RecNetKaon(0),
  Profile_ptSTAR_Term1C2RecNetCharge(0),
  Profile_ptSTAR_Term2C2RecNetCharge(0),
  Profile_ptSTAR_Term1C2RecNetProton(0),
  Profile_ptSTAR_Term2C2RecNetProton(0),
  Profile_ptmax2_CorrectedNetKaon(0),
  Profile_ptmax2_CorrectedNetProton(0),
  Profile_ptmax2_CorrectedNetCharge(0),
  Profile_ptmax2_CorrectedNetChargeNetKaon_term1(0),
  Profile_ptmax2_CorrectedNetChargeNetKaon_term2(0),
  Profile_ptmax2_CorrectedNetChargeNetKaon_term3(0),
  Profile_ptmax2_CorrectedNetChargeNetProton_term1(0),
  Profile_ptmax2_CorrectedNetChargeNetProton_term2(0),
  Profile_ptmax2_CorrectedNetChargeNetProton_term3(0),
  Profile_ptmax2_CorrectedNetKaonNetProton(0),
  Profile_ptmax2_Term1C2CorrectedNetKaon(0),
  Profile_ptmax2_Term2C2CorrectedNetKaon(0),
  Profile_ptmax2_Term3C2CorrectedNetKaon(0),
  Profile_ptmax2_Term4C2CorrectedNetKaon(0),
  Profile_ptmax2_Term1C2CorrectedNetCharge(0),
  Profile_ptmax2_Term2C2CorrectedNetCharge(0),
  Profile_ptmax2_Term3C2CorrectedNetCharge(0),
  Profile_ptmax2_Term4C2CorrectedNetCharge(0),
  Profile_ptmax2_Term1C2CorrectedNetProton(0),
  Profile_ptmax2_Term2C2CorrectedNetProton(0),
  Profile_ptmax2_Term3C2CorrectedNetProton(0),
  Profile_ptmax2_Term4C2CorrectedNetProton(0),
  Profile_ptSTAR_CorrectedNetKaon(0),
  Profile_ptSTAR_CorrectedNetProton(0),
  Profile_ptSTAR_CorrectedNetCharge(0),
  Profile_ptSTAR_CorrectedNetChargeNetKaon_term1(0),
  Profile_ptSTAR_CorrectedNetChargeNetKaon_term2(0),
  Profile_ptSTAR_CorrectedNetChargeNetKaon_term3(0),
  Profile_ptSTAR_CorrectedNetChargeNetProton_term1(0),
  Profile_ptSTAR_CorrectedNetChargeNetProton_term2(0),
  Profile_ptSTAR_CorrectedNetChargeNetProton_term3(0),
  Profile_ptSTAR_CorrectedNetKaonNetProton(0),
  Profile_ptSTAR_Term1C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term2C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term3C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term4C2CorrectedNetKaon(0),
  Profile_ptSTAR_Term1C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term2C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term3C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term4C2CorrectedNetCharge(0),
  Profile_ptSTAR_Term1C2CorrectedNetProton(0),
  Profile_ptSTAR_Term2C2CorrectedNetProton(0),
  Profile_ptSTAR_Term3C2CorrectedNetProton(0),
  Profile_ptSTAR_Term4C2CorrectedNetProton(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaPionCut(0),
  fPIDnSigmaKaonCut(0),
  fPIDnSigmaProtonCut(0),
  fTPCcrossedrows(0),
  fEtaMax(0),
  fPileupCutVal(0),
  fCentralityEstimator_flag(0),
  f2Dhist_nSigmaTPC_pion(0),
  f2Dhist_nSigmaTPC_kaon(0),
  f2Dhist_nSigmaTPC_proton(0),
  f2Dhist_nSigmaTOF_pion(0),
  f2Dhist_nSigmaTOF_kaon(0),
  f2Dhist_nSigmaTOF_proton(0),
  f2Dhist_nSigmaTPCplusTOF_pion(0),
  f2Dhist_nSigmaTPCplusTOF_kaon(0),
  f2Dhist_nSigmaTPCplusTOF_proton(0),
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
  hist_afterCut_TPCncrossedrows(0)
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
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPbPb_v2::~AliAnalysisTaskCorrPbPb_v2()  {

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
  /*
  if (fTreeEvent){
    delete fTreeEvent;
    fTreeEvent = 0x0;
  }
  */
  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }


}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPb_v2::UserCreateOutputObjects()  {
    
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

    
    //QA Plots of Event Selection
    fAODeventCuts.AddQAplotsToList(fQAList,kTRUE);
    fAODeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE, fPileupCutVal);
    
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    // hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    

    //Number of Kaon finally getting selected with specified cuts event wise
    hNumberOfKaonPlus     = new TH1D ("hNumberOfKaonPlus","",3000,0,3000);
    hNumberOfKaonMinus = new TH1D ("hNumberOfKaonMinus","",3000,0,3000);
    // hNumberOfKaonPlus     -> Sumw2();
    // hNumberOfKaonMinus -> Sumw2();
    fOutputList -> Add(hNumberOfKaonPlus);
    fOutputList -> Add(hNumberOfKaonMinus);

    //Number of Pion finally getting selected with specified cuts event wise
    hNumberOfPionPlus     = new TH1D ("hNumberOfPionPlus","",3000,0,3000);
    hNumberOfPionMinus = new TH1D ("hNumberOfPionMinus","",3000,0,3000);
    // hNumberOfPionPlus     -> Sumw2();
    // hNumberOfPionMinus -> Sumw2();
    fOutputList -> Add(hNumberOfPionPlus);
    fOutputList -> Add(hNumberOfPionMinus);

    //Number of Proton finally getting selected with specified cuts event wise
    hNumberOfProtonPlus     = new TH1D ("hNumberOfProtonPlus","",3000,0,3000);
    hNumberOfProtonMinus = new TH1D ("hNumberOfProtonMinus","",3000,0,3000);
    // hNumberOfProtonPlus     -> Sumw2();
    // hNumberOfProtonMinus -> Sumw2();
    fOutputList -> Add(hNumberOfProtonPlus);
    fOutputList -> Add(hNumberOfProtonMinus);


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
    // Profiles for analysis

    //Reconstructed
    
    Profile_ptmax2_RecNetKaon = new TProfile ("Profile_ptmax2_RecNetKaon", "Profile_ptmax2_RecNetKaon", 20, 0, 100);
    Profile_ptmax2_RecNetProton = new TProfile ("Profile_ptmax2_RecNetProton", "Profile_ptmax2_RecNetProton", 20, 0, 100);
    Profile_ptmax2_RecNetCharge = new TProfile ("Profile_ptmax2_RecNetCharge", "Profile_ptmax2_RecNetCharge", 20, 0, 100);
    Profile_ptmax2_RecNetChargeNetKaon = new TProfile ("Profile_ptmax2_RecNetChargeNetKaon", "Profile_ptmax2_RecNetChargeNetKaon", 20, 0, 100);
    Profile_ptmax2_RecNetChargeNetProton = new TProfile ("Profile_ptmax2_RecNetChargeNetProton", "Profile_ptmax2_RecNetChargeNetProton", 20, 0, 100);
    Profile_ptmax2_RecNetKaonNetProton = new TProfile ("Profile_ptmax2_RecNetKaonNetProton", "Profile_ptmax2_RecNetKaonNetProton", 20, 0, 100);
    Profile_ptmax2_Term1C2RecNetKaon = new TProfile ("Profile_ptmax2_Term1C2RecNetKaon", "Profile_ptmax2_Term1C2RecNetKaon", 20, 0, 100);
    Profile_ptmax2_Term1C2RecNetCharge = new TProfile ("Profile_ptmax2_Term1C2RecNetCharge", "Profile_ptmax2_Term1C2RecNetCharge", 20, 0, 100);
    Profile_ptmax2_Term1C2RecNetProton = new TProfile ("Profile_ptmax2_Term1C2RecNetProton", "Profile_ptmax2_Term1C2RecNetProton", 20, 0, 100);
    Profile_ptmax2_Term2C2RecNetKaon = new TProfile ("Profile_ptmax2_Term2C2RecNetKaon", "Profile_ptmax2_Term2C2RecNetKaon", 20, 0, 100);
    Profile_ptmax2_Term2C2RecNetCharge = new TProfile ("Profile_ptmax2_Term2C2RecNetCharge", "Profile_ptmax2_Term2C2RecNetCharge", 20, 0, 100);
    Profile_ptmax2_Term2C2RecNetProton = new TProfile ("Profile_ptmax2_Term2C2RecNetProton", "Profile_ptmax2_Term2C2RecNetProton", 20, 0, 100);

    Profile_ptSTAR_RecNetKaon = new TProfile ("Profile_ptSTAR_RecNetKaon", "Profile_ptSTAR_RecNetKaon", 20, 0, 100);
    Profile_ptSTAR_RecNetProton = new TProfile ("Profile_ptSTAR_RecNetProton", "Profile_ptSTAR_RecNetProton", 20, 0, 100);
    Profile_ptSTAR_RecNetCharge = new TProfile ("Profile_ptSTAR_RecNetCharge", "Profile_ptSTAR_RecNetCharge", 20, 0, 100);
    Profile_ptSTAR_RecNetChargeNetKaon = new TProfile ("Profile_ptSTAR_RecNetChargeNetKaon", "Profile_ptSTAR_RecNetChargeNetKaon", 20, 0, 100);
    Profile_ptSTAR_RecNetChargeNetProton = new TProfile ("Profile_ptSTAR_RecNetChargeNetProton", "Profile_ptSTAR_RecNetChargeNetProton", 20, 0, 100);
    Profile_ptSTAR_RecNetKaonNetProton = new TProfile ("Profile_ptSTAR_RecNetKaonNetProton", "Profile_ptSTAR_RecNetKaonNetProton", 20, 0, 100);
    Profile_ptSTAR_Term1C2RecNetKaon = new TProfile ("Profile_ptSTAR_Term1C2RecNetKaon", "Profile_ptSTAR_Term1C2RecNetKaon", 20, 0, 100);
    Profile_ptSTAR_Term1C2RecNetCharge = new TProfile ("Profile_ptSTAR_Term1C2RecNetCharge", "Profile_ptSTAR_Term1C2RecNetCharge", 20, 0, 100);
    Profile_ptSTAR_Term1C2RecNetProton = new TProfile ("Profile_ptSTAR_Term1C2RecNetProton", "Profile_ptSTAR_Term1C2RecNetProton", 20, 0, 100);
    Profile_ptSTAR_Term2C2RecNetKaon = new TProfile ("Profile_ptSTAR_Term2C2RecNetKaon", "Profile_ptSTAR_Term2C2RecNetKaon", 20, 0, 100);
    Profile_ptSTAR_Term2C2RecNetCharge = new TProfile ("Profile_ptSTAR_Term2C2RecNetCharge", "Profile_ptSTAR_Term2C2RecNetCharge", 20, 0, 100);
    Profile_ptSTAR_Term2C2RecNetProton = new TProfile ("Profile_ptSTAR_Term2C2RecNetProton", "Profile_ptSTAR_Term2C2RecNetProton", 20, 0, 100);

    fTreeList->Add(Profile_ptmax2_RecNetKaon);
    fTreeList->Add(Profile_ptmax2_RecNetProton);
    fTreeList->Add(Profile_ptmax2_RecNetCharge);
    fTreeList->Add(Profile_ptmax2_RecNetChargeNetKaon);
    fTreeList->Add(Profile_ptmax2_RecNetChargeNetProton);
    fTreeList->Add(Profile_ptmax2_RecNetKaonNetProton);
    fTreeList->Add(Profile_ptmax2_Term1C2RecNetKaon);
    fTreeList->Add(Profile_ptmax2_Term2C2RecNetKaon);
    fTreeList->Add(Profile_ptmax2_Term1C2RecNetCharge);
    fTreeList->Add(Profile_ptmax2_Term2C2RecNetCharge);
    fTreeList->Add(Profile_ptmax2_Term1C2RecNetProton);
    fTreeList->Add(Profile_ptmax2_Term2C2RecNetProton);

    fTreeList->Add(Profile_ptSTAR_RecNetKaon);
    fTreeList->Add(Profile_ptSTAR_RecNetProton);
    fTreeList->Add(Profile_ptSTAR_RecNetCharge);
    fTreeList->Add(Profile_ptSTAR_RecNetChargeNetKaon);
    fTreeList->Add(Profile_ptSTAR_RecNetChargeNetProton);
    fTreeList->Add(Profile_ptSTAR_RecNetKaonNetProton);
    fTreeList->Add(Profile_ptSTAR_Term1C2RecNetKaon);
    fTreeList->Add(Profile_ptSTAR_Term2C2RecNetKaon);
    fTreeList->Add(Profile_ptSTAR_Term1C2RecNetCharge);
    fTreeList->Add(Profile_ptSTAR_Term2C2RecNetCharge);
    fTreeList->Add(Profile_ptSTAR_Term1C2RecNetProton);
    fTreeList->Add(Profile_ptSTAR_Term2C2RecNetProton);
   
    
    //Corrected

    //0.2 < pT < 2.0
    
    Profile_ptmax2_CorrectedNetKaon = new TProfile ("Profile_ptmax2_CorrectedNetKaon", "Profile_ptmax2_CorrectedNetKaon", 20, 0, 100);
    Profile_ptmax2_CorrectedNetProton = new TProfile ("Profile_ptmax2_CorrectedNetProton", "Profile_ptmax2_CorrectedNetProton", 20, 0, 100);
    Profile_ptmax2_CorrectedNetCharge = new TProfile ("Profile_ptmax2_CorrectedNetCharge", "Profile_ptmax2_CorrectedNetCharge", 20, 0, 100);

    Profile_ptmax2_CorrectedNetChargeNetKaon_term1 = new TProfile ("Profile_ptmax2_CorrectedNetChargeNetKaon_term1", "Profile_ptmax2_CorrectedNetChargeNetKaon_term1", 20, 0, 100);
    Profile_ptmax2_CorrectedNetChargeNetKaon_term2 = new TProfile ("Profile_ptmax2_CorrectedNetChargeNetKaon_term2", "Profile_ptmax2_CorrectedNetChargeNetKaon_term2", 20, 0, 100);
    Profile_ptmax2_CorrectedNetChargeNetKaon_term3 = new TProfile ("Profile_ptmax2_CorrectedNetChargeNetKaon_term3", "Profile_ptmax2_CorrectedNetChargeNetKaon_term3", 20, 0, 100);

    Profile_ptmax2_CorrectedNetChargeNetProton_term1 = new TProfile ("Profile_ptmax2_CorrectedNetChargeNetProton_term1", "Profile_ptmax2_CorrectedNetChargeNetProton_term1", 20, 0, 100);
    Profile_ptmax2_CorrectedNetChargeNetProton_term2 = new TProfile ("Profile_ptmax2_CorrectedNetChargeNetProton_term2", "Profile_ptmax2_CorrectedNetChargeNetProton_term2", 20, 0, 100);
    Profile_ptmax2_CorrectedNetChargeNetProton_term3 = new TProfile ("Profile_ptmax2_CorrectedNetChargeNetProton_term3", "Profile_ptmax2_CorrectedNetChargeNetProton_term3", 20, 0, 100);

    Profile_ptmax2_CorrectedNetKaonNetProton = new TProfile ("Profile_ptmax2_CorrectedNetKaonNetProton", "Profile_ptmax2_CorrectedNetKaonNetProton", 20, 0, 100);

    Profile_ptmax2_Term1C2CorrectedNetCharge = new TProfile ("Profile_ptmax2_Term1C2CorrectedNetCharge", "Profile_ptmax2_Term1C2CorrectedNetCharge", 20, 0, 100);
    Profile_ptmax2_Term2C2CorrectedNetCharge = new TProfile ("Profile_ptmax2_Term2C2CorrectedNetCharge", "Profile_ptmax2_Term2C2CorrectedNetCharge", 20, 0, 100);
    Profile_ptmax2_Term3C2CorrectedNetCharge = new TProfile ("Profile_ptmax2_Term3C2CorrectedNetCharge", "Profile_ptmax2_Term3C2CorrectedNetCharge", 20, 0, 100);
    Profile_ptmax2_Term4C2CorrectedNetCharge = new TProfile ("Profile_ptmax2_Term4C2CorrectedNetCharge", "Profile_ptmax2_Term4C2CorrectedNetCharge", 20, 0, 100);

    Profile_ptmax2_Term1C2CorrectedNetKaon = new TProfile ("Profile_ptmax2_Term1C2CorrectedNetKaon", "Profile_ptmax2_Term1C2CorrectedNetKaon", 20, 0, 100);
    Profile_ptmax2_Term2C2CorrectedNetKaon = new TProfile ("Profile_ptmax2_Term2C2CorrectedNetKaon", "Profile_ptmax2_Term2C2CorrectedNetKaon", 20, 0, 100);
    Profile_ptmax2_Term3C2CorrectedNetKaon = new TProfile ("Profile_ptmax2_Term3C2CorrectedNetKaon", "Profile_ptmax2_Term3C2CorrectedNetKaon", 20, 0, 100);
    Profile_ptmax2_Term4C2CorrectedNetKaon = new TProfile ("Profile_ptmax2_Term4C2CorrectedNetKaon", "Profile_ptmax2_Term4C2CorrectedNetKaon", 20, 0, 100);

    Profile_ptmax2_Term1C2CorrectedNetProton = new TProfile ("Profile_ptmax2_Term1C2CorrectedNetProton", "Profile_ptmax2_Term1C2CorrectedNetProton", 20, 0, 100);
    Profile_ptmax2_Term2C2CorrectedNetProton = new TProfile ("Profile_ptmax2_Term2C2CorrectedNetProton", "Profile_ptmax2_Term2C2CorrectedNetProton", 20, 0, 100);
    Profile_ptmax2_Term3C2CorrectedNetProton = new TProfile ("Profile_ptmax2_Term3C2CorrectedNetProton", "Profile_ptmax2_Term3C2CorrectedNetProton", 20, 0, 100);
    Profile_ptmax2_Term4C2CorrectedNetProton = new TProfile ("Profile_ptmax2_Term4C2CorrectedNetProton", "Profile_ptmax2_Term4C2CorrectedNetProton", 20, 0, 100);

    fTreeList->Add(Profile_ptmax2_CorrectedNetKaon);
    fTreeList->Add(Profile_ptmax2_CorrectedNetProton);
    fTreeList->Add(Profile_ptmax2_CorrectedNetCharge);
    fTreeList->Add(Profile_ptmax2_CorrectedNetChargeNetKaon_term1);
    fTreeList->Add(Profile_ptmax2_CorrectedNetChargeNetKaon_term2);
    fTreeList->Add(Profile_ptmax2_CorrectedNetChargeNetKaon_term3);
    fTreeList->Add(Profile_ptmax2_CorrectedNetChargeNetProton_term1);
    fTreeList->Add(Profile_ptmax2_CorrectedNetChargeNetProton_term2);
    fTreeList->Add(Profile_ptmax2_CorrectedNetChargeNetProton_term3);
    fTreeList->Add(Profile_ptmax2_CorrectedNetKaonNetProton);
    fTreeList->Add(Profile_ptmax2_Term1C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptmax2_Term2C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptmax2_Term3C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptmax2_Term4C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptmax2_Term1C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptmax2_Term2C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptmax2_Term3C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptmax2_Term4C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptmax2_Term1C2CorrectedNetProton);
    fTreeList->Add(Profile_ptmax2_Term2C2CorrectedNetProton);
    fTreeList->Add(Profile_ptmax2_Term3C2CorrectedNetProton);
    fTreeList->Add(Profile_ptmax2_Term4C2CorrectedNetProton);


    //pt STAR cut
    
    Profile_ptSTAR_CorrectedNetKaon = new TProfile ("Profile_ptSTAR_CorrectedNetKaon", "Profile_ptSTAR_CorrectedNetKaon", 20, 0, 100);
    Profile_ptSTAR_CorrectedNetProton = new TProfile ("Profile_ptSTAR_CorrectedNetProton", "Profile_ptSTAR_CorrectedNetProton", 20, 0, 100);
    Profile_ptSTAR_CorrectedNetCharge = new TProfile ("Profile_ptSTAR_CorrectedNetCharge", "Profile_ptSTAR_CorrectedNetCharge", 20, 0, 100);

    Profile_ptSTAR_CorrectedNetChargeNetKaon_term1 = new TProfile ("Profile_ptSTAR_CorrectedNetChargeNetKaon_term1", "Profile_ptSTAR_CorrectedNetChargeNetKaon_term1", 20, 0, 100);
    Profile_ptSTAR_CorrectedNetChargeNetKaon_term2 = new TProfile ("Profile_ptSTAR_CorrectedNetChargeNetKaon_term2", "Profile_ptSTAR_CorrectedNetChargeNetKaon_term2", 20, 0, 100);
    Profile_ptSTAR_CorrectedNetChargeNetKaon_term3 = new TProfile ("Profile_ptSTAR_CorrectedNetChargeNetKaon_term3", "Profile_ptSTAR_CorrectedNetChargeNetKaon_term3", 20, 0, 100);

    Profile_ptSTAR_CorrectedNetChargeNetProton_term1 = new TProfile ("Profile_ptSTAR_CorrectedNetChargeNetProton_term1", "Profile_ptSTAR_CorrectedNetChargeNetProton_term1", 20, 0, 100);
    Profile_ptSTAR_CorrectedNetChargeNetProton_term2 = new TProfile ("Profile_ptSTAR_CorrectedNetChargeNetProton_term2", "Profile_ptSTAR_CorrectedNetChargeNetProton_term2", 20, 0, 100);
    Profile_ptSTAR_CorrectedNetChargeNetProton_term3 = new TProfile ("Profile_ptSTAR_CorrectedNetChargeNetProton_term3", "Profile_ptSTAR_CorrectedNetChargeNetProton_term3", 20, 0, 100);

    Profile_ptSTAR_CorrectedNetKaonNetProton = new TProfile ("Profile_ptSTAR_CorrectedNetKaonNetProton", "Profile_ptSTAR_CorrectedNetKaonNetProton", 20, 0, 100);

    Profile_ptSTAR_Term1C2CorrectedNetCharge = new TProfile ("Profile_ptSTAR_Term1C2CorrectedNetCharge", "Profile_ptSTAR_Term1C2CorrectedNetCharge", 20, 0, 100);
    Profile_ptSTAR_Term2C2CorrectedNetCharge = new TProfile ("Profile_ptSTAR_Term2C2CorrectedNetCharge", "Profile_ptSTAR_Term2C2CorrectedNetCharge", 20, 0, 100);
    Profile_ptSTAR_Term3C2CorrectedNetCharge = new TProfile ("Profile_ptSTAR_Term3C2CorrectedNetCharge", "Profile_ptSTAR_Term3C2CorrectedNetCharge", 20, 0, 100);
    Profile_ptSTAR_Term4C2CorrectedNetCharge = new TProfile ("Profile_ptSTAR_Term4C2CorrectedNetCharge", "Profile_ptSTAR_Term4C2CorrectedNetCharge", 20, 0, 100);

    Profile_ptSTAR_Term1C2CorrectedNetKaon = new TProfile ("Profile_ptSTAR_Term1C2CorrectedNetKaon", "Profile_ptSTAR_Term1C2CorrectedNetKaon", 20, 0, 100);
    Profile_ptSTAR_Term2C2CorrectedNetKaon = new TProfile ("Profile_ptSTAR_Term2C2CorrectedNetKaon", "Profile_ptSTAR_Term2C2CorrectedNetKaon", 20, 0, 100);
    Profile_ptSTAR_Term3C2CorrectedNetKaon = new TProfile ("Profile_ptSTAR_Term3C2CorrectedNetKaon", "Profile_ptSTAR_Term3C2CorrectedNetKaon", 20, 0, 100);
    Profile_ptSTAR_Term4C2CorrectedNetKaon = new TProfile ("Profile_ptSTAR_Term4C2CorrectedNetKaon", "Profile_ptSTAR_Term4C2CorrectedNetKaon", 20, 0, 100);

    Profile_ptSTAR_Term1C2CorrectedNetProton = new TProfile ("Profile_ptSTAR_Term1C2CorrectedNetProton", "Profile_ptSTAR_Term1C2CorrectedNetProton", 20, 0, 100);
    Profile_ptSTAR_Term2C2CorrectedNetProton = new TProfile ("Profile_ptSTAR_Term2C2CorrectedNetProton", "Profile_ptSTAR_Term2C2CorrectedNetProton", 20, 0, 100);
    Profile_ptSTAR_Term3C2CorrectedNetProton = new TProfile ("Profile_ptSTAR_Term3C2CorrectedNetProton", "Profile_ptSTAR_Term3C2CorrectedNetProton", 20, 0, 100);
    Profile_ptSTAR_Term4C2CorrectedNetProton = new TProfile ("Profile_ptSTAR_Term4C2CorrectedNetProton", "Profile_ptSTAR_Term4C2CorrectedNetProton", 20, 0, 100);

    fTreeList->Add(Profile_ptSTAR_CorrectedNetKaon);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetProton);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetCharge);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetChargeNetKaon_term1);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetChargeNetKaon_term2);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetChargeNetKaon_term3);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetChargeNetProton_term1);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetChargeNetProton_term2);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetChargeNetProton_term3);
    fTreeList->Add(Profile_ptSTAR_CorrectedNetKaonNetProton);
    fTreeList->Add(Profile_ptSTAR_Term1C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptSTAR_Term2C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptSTAR_Term3C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptSTAR_Term4C2CorrectedNetKaon);
    fTreeList->Add(Profile_ptSTAR_Term1C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptSTAR_Term2C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptSTAR_Term3C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptSTAR_Term4C2CorrectedNetCharge);
    fTreeList->Add(Profile_ptSTAR_Term1C2CorrectedNetProton);
    fTreeList->Add(Profile_ptSTAR_Term2C2CorrectedNetProton);
    fTreeList->Add(Profile_ptSTAR_Term3C2CorrectedNetProton);
    fTreeList->Add(Profile_ptSTAR_Term4C2CorrectedNetProton);


    
    
    //TTree object to store variables
    /*
    fTreeEvent = new TTree("fTreeEvent","Event Tree");
    fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    //reconstructed
    fTreeEvent->Branch("fNoKaonPlus_ptmax2", &fNoKaonPlus_ptmax2, "fNoKaonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoKaonMinus_ptmax2", &fNoKaonMinus_ptmax2, "fNoKaonMinus_ptmax2/F");
    fTreeEvent->Branch("fNoPionPlus_ptmax2", &fNoPionPlus_ptmax2, "fNoPionPlus_ptmax2/F");
    fTreeEvent->Branch("fNoPionMinus_ptmax2", &fNoPionMinus_ptmax2, "fNoPionMinus_ptmax2/F");
    fTreeEvent->Branch("fNoProtonPlus_ptmax2", &fNoProtonPlus_ptmax2, "fNoProtonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoProtonMinus_ptmax2", &fNoProtonMinus_ptmax2, "fNoProtonMinus_ptmax2/F");

    fTreeEvent->Branch("fNoKaonPlus_ptmax3", &fNoKaonPlus_ptmax3, "fNoKaonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoKaonMinus_ptmax3", &fNoKaonMinus_ptmax3, "fNoKaonMinus_ptmax3/F");
    fTreeEvent->Branch("fNoPionPlus_ptmax3", &fNoPionPlus_ptmax3, "fNoPionPlus_ptmax3/F");
    fTreeEvent->Branch("fNoPionMinus_ptmax3", &fNoPionMinus_ptmax3, "fNoPionMinus_ptmax3/F");
    fTreeEvent->Branch("fNoProtonPlus_ptmax3", &fNoProtonPlus_ptmax3, "fNoProtonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoProtonMinus_ptmax3", &fNoProtonMinus_ptmax3, "fNoProtonMinus_ptmax3/F");

    //Reconstred tracks corrected
    fTreeEvent->Branch("fCorrectedNoKaonPlus_ptmax2", &fCorrectedNoKaonPlus_ptmax2, "fCorrectedNoKaonPlus_ptmax2/F");
    fTreeEvent->Branch("fCorrectedNoKaonMinus_ptmax2", &fCorrectedNoKaonMinus_ptmax2, "fCorrectedNoKaonMinus_ptmax2/F");
    fTreeEvent->Branch("fCorrectedNoPionPlus_ptmax2", &fCorrectedNoPionPlus_ptmax2, "fCorrectedNoPionPlus_ptmax2/F");
    fTreeEvent->Branch("fCorrectedNoPionMinus_ptmax2", &fCorrectedNoPionMinus_ptmax2, "fCorrectedNoPionMinus_ptmax2/F");
    fTreeEvent->Branch("fCorrectedNoProtonPlus_ptmax2", &fCorrectedNoProtonPlus_ptmax2, "fCorrectedNoProtonPlus_ptmax2/F");
    fTreeEvent->Branch("fCorrectedNoProtonMinus_ptmax2", &fCorrectedNoProtonMinus_ptmax2, "fCorrectedNoProtonMinus_ptmax2/F");

    fTreeEvent->Branch("fCorrectedNoKaonPlus_ptmax3", &fCorrectedNoKaonPlus_ptmax3, "fCorrectedNoKaonPlus_ptmax3/F");
    fTreeEvent->Branch("fCorrectedNoKaonMinus_ptmax3", &fCorrectedNoKaonMinus_ptmax3, "fCorrectedNoKaonMinus_ptmax3/F");
    fTreeEvent->Branch("fCorrectedNoPionPlus_ptmax3", &fCorrectedNoPionPlus_ptmax3, "fCorrectedNoPionPlus_ptmax3/F");
    fTreeEvent->Branch("fCorrectedNoPionMinus_ptmax3", &fCorrectedNoPionMinus_ptmax3, "fCorrectedNoPionMinus_ptmax3/F");
    fTreeEvent->Branch("fCorrectedNoProtonPlus_ptmax3", &fCorrectedNoProtonPlus_ptmax3, "fCorrectedNoProtonPlus_ptmax3/F");
    fTreeEvent->Branch("fCorrectedNoProtonMinus_ptmax3", &fCorrectedNoProtonMinus_ptmax3, "fCorrectedNoProtonMinus_ptmax3/F");

    fTreeEvent->Branch("fEffSqrFactrKaonPlus_ptmax2", &fEffSqrFactrKaonPlus_ptmax2, "fEffSqrFactrKaonPlus_ptmax2/F");
    fTreeEvent->Branch("fEffSqrFactrKaonMinus_ptmax2", &fEffSqrFactrKaonMinus_ptmax2, "fEffSqrFactrKaonMinus_ptmax2/F");
    fTreeEvent->Branch("fEffSqrFactrPionPlus_ptmax2", &fEffSqrFactrPionPlus_ptmax2, "fEffSqrFactrPionPlus_ptmax2/F");
    fTreeEvent->Branch("fEffSqrFactrPionMinus_ptmax2", &fEffSqrFactrPionMinus_ptmax2, "fEffSqrFactrPionMinus_ptmax2/F");
    fTreeEvent->Branch("fEffSqrFactrProtonPlus_ptmax2", &fEffSqrFactrProtonPlus_ptmax2, "fEffSqrFactrProtonPlus_ptmax2/F");
    fTreeEvent->Branch("fEffSqrFactrProtonMinus_ptmax2", &fEffSqrFactrProtonMinus_ptmax2, "fEffSqrFactrProtonMinus_ptmax2/F");
    
    fTreeEvent->Branch("fEffSqrFactrKaonPlus_ptmax3", &fEffSqrFactrKaonPlus_ptmax3, "fEffSqrFactrKaonPlus_ptmax3/F");
    fTreeEvent->Branch("fEffSqrFactrKaonMinus_ptmax3", &fEffSqrFactrKaonMinus_ptmax3, "fEffSqrFactrKaonMinus_ptmax3/F");
    fTreeEvent->Branch("fEffSqrFactrPionPlus_ptmax3", &fEffSqrFactrPionPlus_ptmax3, "fEffSqrFactrPionPlus_ptmax3/F");
    fTreeEvent->Branch("fEffSqrFactrPionMinus_ptmax3", &fEffSqrFactrPionMinus_ptmax3, "fEffSqrFactrPionMinus_ptmax3/F");
    fTreeEvent->Branch("fEffSqrFactrProtonPlus_ptmax3", &fEffSqrFactrProtonPlus_ptmax3, "fEffSqrFactrProtonPlus_ptmax3/F");
    fTreeEvent->Branch("fEffSqrFactrProtonMinus_ptmax3", &fEffSqrFactrProtonMinus_ptmax3, "fEffSqrFactrProtonMinus_ptmax3/F");

    fTreeList->Add(fTreeEvent);
    */
    
    /*
    
    //Track Cuts Objects
    if(!fESDtrackCuts )
      {
	fESDtrackCuts = new AliESDtrackCuts();
	fESDtrackCuts->SetMinNClustersTPC(70);
	fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
	fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCuts->SetRequireTPCRefit(kTRUE);
	fESDtrackCuts->SetPtRange(0.15,1e10);
	fESDtrackCuts->SetEtaRange(-0.8,0.8);
      }
    if(!fESDtrackCuts_primary)
      {
	fESDtrackCuts_primary = new AliESDtrackCuts();
	fESDtrackCuts_primary = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);
	fESDtrackCuts_primary->SetPtRange(0.3,5.0);
	fESDtrackCuts_primary->SetEtaRange(-0.8,+0.8);
      }
    */
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeList);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPb_v2::UserExec(Option_t *)  {
  
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

    
    //"fTreeEvent" Tree Variable
    //fTreeVariableCentrality = lV0M;

    
    
    //Initialize number 0f K+, K-, pi+, pi-, p and p-bar per event

    //raw
    Float_t no_KaonPlus_perevent = 0;
    Float_t no_KaonMinus_perevent = 0;
    Float_t no_ProtonPlus_perevent = 0;
    Float_t no_ProtonMinus_perevent = 0;
    Float_t no_PionPlus_perevent = 0;
    Float_t no_PionMinus_perevent = 0;

    Float_t no_KaonPlus_perevent_ptmax2 = 0;
    Float_t no_KaonMinus_perevent_ptmax2 = 0;
    Float_t no_ProtonPlus_perevent_ptmax2 = 0;
    Float_t no_ProtonMinus_perevent_ptmax2 = 0;
    Float_t no_PionPlus_perevent_ptmax2 = 0;
    Float_t no_PionMinus_perevent_ptmax2 = 0;

    //efficiency corrected
    Float_t no_KaonPlus_perevent_corrected = 0;
    Float_t no_KaonMinus_perevent_corrected = 0;
    Float_t no_ProtonPlus_perevent_corrected = 0;
    Float_t no_ProtonMinus_perevent_corrected = 0;
    Float_t no_PionPlus_perevent_corrected = 0;
    Float_t no_PionMinus_perevent_corrected = 0;

    Float_t no_KaonPlus_perevent_ptmax2_corrected = 0;
    Float_t no_KaonMinus_perevent_ptmax2_corrected = 0;
    Float_t no_ProtonPlus_perevent_ptmax2_corrected = 0;
    Float_t no_ProtonMinus_perevent_ptmax2_corrected = 0;
    Float_t no_PionPlus_perevent_ptmax2_corrected = 0;
    Float_t no_PionMinus_perevent_ptmax2_corrected = 0;

    Float_t noByEffSquare_KaonPlus_perevent_corrected = 0;
    Float_t noByEffSquare_KaonMinus_perevent_corrected = 0;
    Float_t noByEffSquare_ProtonPlus_perevent_corrected = 0;
    Float_t noByEffSquare_ProtonMinus_perevent_corrected = 0;
    Float_t noByEffSquare_PionPlus_perevent_corrected = 0;
    Float_t noByEffSquare_PionMinus_perevent_corrected = 0;

    Float_t noByEffSquare_KaonPlus_perevent_ptmax2_corrected = 0;
    Float_t noByEffSquare_KaonMinus_perevent_ptmax2_corrected = 0;
    Float_t noByEffSquare_ProtonPlus_perevent_ptmax2_corrected = 0;
    Float_t noByEffSquare_ProtonMinus_perevent_ptmax2_corrected = 0;
    Float_t noByEffSquare_PionPlus_perevent_ptmax2_corrected = 0;
    Float_t noByEffSquare_PionMinus_perevent_ptmax2_corrected = 0;


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Int_t ptBinNo = 0;
    Double_t BinCont = 0;
    Double_t EffWgt = 0;

    if(fListTRKCorr) 
      {
	cout<<"## Got efficiency histograms..## "<<endl;
	GetMCEffCorrectionHist();
      }
    else
      cout<<"################ No histograms found ############### "<<endl;

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




    //Loop on reconstructed tracks
    
    for(Int_t itr=0; itr < fAODevent->GetNumberOfTracks(); itr++)
      {
	
	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;

	//getting all track variables
	
	Double_t trkPt = aodtrack->Pt();
	Double_t trkPhi = aodtrack->Phi();
	Double_t trkEta = aodtrack->Eta();
	Double_t trkCharge = aodtrack->Charge();
	Float_t trkDCAxy;
	Float_t trkDCAz;
	aodtrack->GetImpactParameters(trkDCAxy, trkDCAz);
	Double_t trkTPCNCls = aodtrack->GetTPCNcls();
	Double_t trkChi2PerNDF = aodtrack->Chi2perNDF();
	Double_t trkITSchi2 = aodtrack->GetITSchi2();
	Int_t trkITSNcls = aodtrack->GetITSNcls();
	Double_t trkITSchi2perNcls = trkITSchi2/trkITSNcls;
	Double_t trkTPCchi2perNcls = aodtrack->GetTPCchi2perCluster();
	Double_t trkTPCcrossedrows = aodtrack->GetTPCCrossedRows();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Histograms filled befor applying track cut
	hist_beforeCut_DCAxy->Fill(trkDCAxy);
	hist_beforeCut_DCAz->Fill(trkDCAz);
	hist_beforeCut_eta->Fill(trkEta);
	hist_beforeCut_chi2perTPCclstr->Fill(trkTPCchi2perNcls);
	hist_beforeCut_chi2perITSclstr->Fill(trkITSchi2perNcls);
	hist_beforeCut_TPCncrossedrows->Fill(trkTPCcrossedrows);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	//All track cuts
	
	//Track selectionL FilterBit : default is FB96
	if(!aodtrack->TestFilterBit(fFBNo))  continue;

	//cuts on TPCchi2perClstr and ITSchi2perClstr and TPC nCrossedRows
	if (trkTPCcrossedrows < fTPCcrossedrows) continue;
	if (trkTPCchi2perNcls > fChi2TPC) continue;
	if (trkITSchi2perNcls > fChi2ITS) continue;


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//PID nSigma histograms

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
	f2Dhist_nSigmaTPC_pion->Fill(fTPCnSigma_Pion, trkPt);
	f2Dhist_nSigmaTPC_kaon->Fill(fTPCnSigma_Kaon, trkPt);
	f2Dhist_nSigmaTPC_proton->Fill(fTPCnSigma_Proton, trkPt);
	//TOF
	f2Dhist_nSigmaTOF_pion->Fill(fTOFnSigma_Pion, trkPt);
	f2Dhist_nSigmaTOF_kaon->Fill(fTOFnSigma_Kaon, trkPt);
	f2Dhist_nSigmaTOF_proton->Fill(fTOFnSigma_Proton, trkPt);

	Int_t flaggg = 0;
	if(trkPt >= 0.5 )
	  {
	    //rejection
	    
	    if (TMath::Abs(fTPCplusTOFnSigma_Pion) < 3.0)
	      flaggg += 1;
	    if (TMath::Abs(fTPCplusTOFnSigma_Proton) < 3.0)
	      flaggg += 1;
	    if (TMath::Abs(fTPCplusTOFnSigma_Kaon) < 3.0)
	      flaggg += 1;
	  }
	if (flaggg > 1) continue;

	//TPC+TOF
	f2Dhist_nSigmaTPCplusTOF_pion->Fill(fTPCplusTOFnSigma_Pion, trkPt);
	f2Dhist_nSigmaTPCplusTOF_kaon->Fill(fTPCplusTOFnSigma_Kaon, trkPt);
	f2Dhist_nSigmaTPCplusTOF_proton->Fill(fTPCplusTOFnSigma_Proton, trkPt);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



	//Kinematic cuts on pT and Eta
	if (TMath::Abs(trkEta) > fEtaMax) continue;
	if (trkPt < 0.2) continue;
	if (trkPt > 3.0) continue;


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Histograms filled after applying track cut
	hist_afterCut_DCAxy->Fill(trkDCAxy);
	hist_afterCut_DCAz->Fill(trkDCAz);
	hist_afterCut_eta->Fill(trkEta);
	hist_afterCut_chi2perTPCclstr->Fill(trkTPCchi2perNcls);
	hist_afterCut_chi2perITSclstr->Fill(trkITSchi2perNcls);
	hist_afterCut_TPCncrossedrows->Fill(trkTPCcrossedrows);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



	//PID selection
	Bool_t IsKaon = KaonSelector(track, fPIDnSigmaKaonCut);
	Bool_t IsPion = PionSelector(track, fPIDnSigmaPionCut);
	Bool_t IsProton = ProtonSelector(track, fPIDnSigmaProtonCut);
	if (!IsKaon && !IsPion && !IsProton) continue;

	//Check if a particle is selected as more than one type
	Int_t flag = 0;
	if(IsKaon) flag+=1;
	if(IsPion) flag+=1;
	if(IsProton) flag+=1;
	if(flag>1)
	  {
	    cout<<"Particle identified as more than on PID: flag= "<<flag<<endl;
	    continue;
	  }
	
	if (trkCharge > 0 && IsKaon)   //K+
	  {
	    if (fEffKaonPlus[centrality_bin])
	      {
		ptBinNo = fEffKaonPlus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffKaonPlus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont;

		if(trkPt > 0.4 && trkPt < 1.6)
		  {
		    no_KaonPlus_perevent += 1.0;
		    no_KaonPlus_perevent_corrected += EffWgt;
		    noByEffSquare_KaonPlus_perevent_corrected += TMath::Power(EffWgt, 2.0);
		  }
	    
		if(trkPt < 2.0)
		  {
		    no_KaonPlus_perevent_ptmax2 += 1.0;
		    no_KaonPlus_perevent_ptmax2_corrected += EffWgt;
		    noByEffSquare_KaonPlus_perevent_ptmax2_corrected += TMath::Power(EffWgt, 2.0);
		  }
	      }
	  }
	  

	if (trkCharge < 0 && IsKaon)   //K-
	  {
	    if (fEffKaonMinus[centrality_bin])
	      {
		ptBinNo = fEffKaonMinus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffKaonMinus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont;

		if(trkPt > 0.4 && trkPt < 1.6)
		  {
		    no_KaonMinus_perevent += 1.0;
		    no_KaonMinus_perevent_corrected += EffWgt;
		    noByEffSquare_KaonMinus_perevent_corrected += TMath::Power(EffWgt, 2.0);
		  }
	   
		if(trkPt < 2.0)
		  {
		    no_KaonMinus_perevent_ptmax2 += 1.0;
		    no_KaonMinus_perevent_ptmax2_corrected += EffWgt;
		    noByEffSquare_KaonMinus_perevent_ptmax2_corrected += TMath::Power(EffWgt, 2.0);
		  }
	      }
	  }

	if (trkCharge > 0 && IsProton)   //proton
	  {
	    if (trkPt > 0.4) // cut for removing protons coming from beam pipe
	      {
		if (fEffProtonPlus[centrality_bin])
		  {
		    ptBinNo = fEffProtonPlus[centrality_bin]->FindBin(trkPt);
		    BinCont = fEffProtonPlus[centrality_bin]->GetBinContent(ptBinNo);
		    if(BinCont!=0) EffWgt = 1.0/BinCont;

		    if(trkPt > 0.4 && trkPt < 1.6)
		      {
			no_ProtonPlus_perevent += 1.0;
			no_ProtonPlus_perevent_corrected += EffWgt;
			noByEffSquare_ProtonPlus_perevent_corrected += TMath::Power(EffWgt, 2.0);
		      }
		    
		    if(trkPt < 2.0)
		      {
			no_ProtonPlus_perevent_ptmax2 += 1.0;
		      	no_ProtonPlus_perevent_ptmax2_corrected += EffWgt;
			noByEffSquare_ProtonPlus_perevent_ptmax2_corrected += TMath::Power(EffWgt, 2.0);
		      }
		  }
	      }
	  }
	  

	if (trkCharge < 0 && IsProton)   //anti-proton
	  {
	    if (trkPt > 0.4) // cut for removing protons coming from beam pipe
	      {
		if (fEffProtonMinus[centrality_bin])
		  {
		    ptBinNo = fEffProtonMinus[centrality_bin]->FindBin(trkPt);
		    BinCont = fEffProtonMinus[centrality_bin]->GetBinContent(ptBinNo);
		    if(BinCont!=0) EffWgt = 1.0/BinCont;

		    if(trkPt > 0.4 && trkPt < 1.6)
		      {
			no_ProtonMinus_perevent += 1.0;
			no_ProtonMinus_perevent_corrected += EffWgt;
			noByEffSquare_ProtonMinus_perevent_corrected += TMath::Power(EffWgt, 2.0);
		      }

		    if(trkPt < 2.0)
		      {
			no_ProtonMinus_perevent_ptmax2 += 1.0;
			no_ProtonMinus_perevent_ptmax2_corrected += EffWgt;
			noByEffSquare_ProtonMinus_perevent_ptmax2_corrected += TMath::Power(EffWgt, 2.0);
		      }
		  }
	      }
	  }

	if (trkCharge > 0 && IsPion)   //pi+
	  {
	    if (fEffPionPlus[centrality_bin])
	      {
		ptBinNo = fEffPionPlus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffPionPlus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont;

		
		if(trkPt > 0.4 && trkPt < 1.6)
		  {
		    no_PionPlus_perevent += 1.0;
		    no_PionPlus_perevent_corrected += EffWgt;
		    noByEffSquare_PionPlus_perevent_corrected += TMath::Power(EffWgt, 2.0);
		  }

		if(trkPt < 2.0)
		  {
		    no_PionPlus_perevent_ptmax2 += 1.0;
		    no_PionPlus_perevent_ptmax2_corrected += EffWgt;
		    noByEffSquare_PionPlus_perevent_ptmax2_corrected += TMath::Power(EffWgt, 2.0);
		  }
	      }
	  }
	  

	if (trkCharge < 0 && IsPion)   //pi-
	  {
	    if (fEffPionMinus[centrality_bin])
	      {
		ptBinNo = fEffPionMinus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffPionMinus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont;
	      
		if(trkPt > 0.4 && trkPt < 1.6)
		  {
		    no_PionMinus_perevent += 1.0;
		    no_PionMinus_perevent_corrected += EffWgt;
		    noByEffSquare_PionMinus_perevent_corrected += TMath::Power(EffWgt, 2.0);
		  }

	
		if(trkPt < 2.0)
		  {
		    no_PionMinus_perevent_ptmax2 += 1.0;
		    no_PionMinus_perevent_ptmax2_corrected += EffWgt;
		    noByEffSquare_PionMinus_perevent_ptmax2_corrected += TMath::Power(EffWgt, 2.0);
		  }
	      }
	  }
      }      
    //end reconstructed track loop
    

    
    //Kaon
    hNumberOfKaonPlus->Fill(no_KaonPlus_perevent_ptmax2);
    hNumberOfKaonMinus->Fill(no_KaonMinus_perevent_ptmax2);
    //Pion
    hNumberOfPionPlus->Fill(no_PionPlus_perevent_ptmax2);
    hNumberOfPionMinus->Fill(no_PionMinus_perevent_ptmax2);
    //Proton
    hNumberOfProtonPlus->Fill(no_ProtonPlus_perevent_ptmax2);
    hNumberOfProtonMinus->Fill(no_ProtonMinus_perevent_ptmax2);

    //++++++++++++++++++++++++++++++++++++++++++++
    //Reconstructed

    /*
    //Kaon
    fNoKaonPlus_ptmax2 = no_KaonPlus_perevent_ptmax2;
    fNoKaonMinus_ptmax2 = no_KaonMinus_perevent_ptmax2;
    fNoKaonPlus_ptmax3 = no_KaonPlus_perevent;
    fNoKaonMinus_ptmax3 = no_KaonMinus_perevent;
    //Pion
    fNoPionPlus_ptmax2 = no_PionPlus_perevent_ptmax2;
    fNoPionMinus_ptmax2 = no_PionMinus_perevent_ptmax2;
    fNoPionPlus_ptmax3 = no_PionPlus_perevent;
    fNoPionMinus_ptmax3 = no_PionMinus_perevent;
    //Proton
    fNoProtonPlus_ptmax2 = no_ProtonPlus_perevent_ptmax2;
    fNoProtonMinus_ptmax2 = no_ProtonMinus_perevent_ptmax2;
    fNoProtonPlus_ptmax3 = no_ProtonPlus_perevent;
    fNoProtonMinus_ptmax3 = no_ProtonMinus_perevent;
    */

    Double_t netKaon_rec_ptmax2 = no_KaonPlus_perevent_ptmax2 - no_KaonMinus_perevent_ptmax2;
    Double_t netCharge_rec_ptmax2 = no_PionPlus_perevent_ptmax2 - no_PionMinus_perevent_ptmax2 + no_KaonPlus_perevent_ptmax2 - no_KaonMinus_perevent_ptmax2 + no_ProtonPlus_perevent_ptmax2 - no_ProtonMinus_perevent_ptmax2;
    Double_t netProton_rec_ptmax2 = no_ProtonPlus_perevent_ptmax2 - no_ProtonMinus_perevent_ptmax2;

    Double_t netKaon_rec_ptSTAR = no_KaonPlus_perevent - no_KaonMinus_perevent;
    Double_t netCharge_rec_ptSTAR = no_PionPlus_perevent - no_PionMinus_perevent + no_KaonPlus_perevent - no_KaonMinus_perevent + no_ProtonPlus_perevent - no_ProtonMinus_perevent;
    Double_t netProton_rec_ptSTAR = no_ProtonPlus_perevent - no_ProtonMinus_perevent;

    Double_t term1_netKaon_rec_ptmax2 = TMath::Power(netKaon_rec_ptmax2, 2.0);
    Double_t term1_netCharge_rec_ptmax2 = TMath::Power(netCharge_rec_ptmax2, 2.0);
    Double_t term1_netProton_rec_ptmax2 = TMath::Power(netProton_rec_ptmax2, 2.0);
    Double_t term2_netKaon_rec_ptmax2 = TMath::Power(netKaon_rec_ptmax2, 1.0);
    Double_t term2_netCharge_rec_ptmax2 = TMath::Power(netCharge_rec_ptmax2, 1.0);
    Double_t term2_netProton_rec_ptmax2 = TMath::Power(netProton_rec_ptmax2, 1.0);

    Double_t term1_netKaon_rec_ptSTAR = TMath::Power(netKaon_rec_ptSTAR, 2.0);
    Double_t term1_netCharge_rec_ptSTAR = TMath::Power(netCharge_rec_ptSTAR, 2.0);
    Double_t term1_netProton_rec_ptSTAR = TMath::Power(netProton_rec_ptSTAR, 2.0);
    Double_t term2_netKaon_rec_ptSTAR = TMath::Power(netKaon_rec_ptSTAR, 1.0);
    Double_t term2_netCharge_rec_ptSTAR = TMath::Power(netCharge_rec_ptSTAR, 1.0);
    Double_t term2_netProton_rec_ptSTAR = TMath::Power(netProton_rec_ptSTAR, 1.0);

    
    
    Profile_ptmax2_RecNetKaon->Fill(lV0M,netKaon_rec_ptmax2);
    Profile_ptmax2_RecNetCharge->Fill(lV0M,netCharge_rec_ptmax2);
    Profile_ptmax2_RecNetProton->Fill(lV0M,netProton_rec_ptmax2);
    Profile_ptmax2_RecNetChargeNetKaon->Fill(lV0M,netCharge_rec_ptmax2*netKaon_rec_ptmax2);
    Profile_ptmax2_RecNetChargeNetProton->Fill(lV0M,netCharge_rec_ptmax2*netProton_rec_ptmax2);
    Profile_ptmax2_RecNetKaonNetProton->Fill(lV0M,netKaon_rec_ptmax2*netProton_rec_ptmax2);

    Profile_ptSTAR_RecNetKaon->Fill(lV0M,netKaon_rec_ptSTAR);
    Profile_ptSTAR_RecNetCharge->Fill(lV0M,netCharge_rec_ptSTAR);
    Profile_ptSTAR_RecNetProton->Fill(lV0M,netProton_rec_ptSTAR);
    Profile_ptSTAR_RecNetChargeNetKaon->Fill(lV0M,netCharge_rec_ptSTAR*netKaon_rec_ptSTAR);
    Profile_ptSTAR_RecNetChargeNetProton->Fill(lV0M,netCharge_rec_ptSTAR*netProton_rec_ptSTAR);
    Profile_ptSTAR_RecNetKaonNetProton->Fill(lV0M,netKaon_rec_ptSTAR*netProton_rec_ptSTAR);

    Profile_ptmax2_Term1C2RecNetKaon->Fill(lV0M, term1_netKaon_rec_ptmax2);
    Profile_ptmax2_Term1C2RecNetCharge->Fill(lV0M, term1_netCharge_rec_ptmax2);
    Profile_ptmax2_Term1C2RecNetProton->Fill(lV0M, term1_netProton_rec_ptmax2);
    Profile_ptmax2_Term2C2RecNetKaon->Fill(lV0M, term2_netKaon_rec_ptmax2);
    Profile_ptmax2_Term2C2RecNetCharge->Fill(lV0M, term2_netCharge_rec_ptmax2);
    Profile_ptmax2_Term2C2RecNetProton->Fill(lV0M, term2_netProton_rec_ptmax2);

    Profile_ptSTAR_Term1C2RecNetKaon->Fill(lV0M, term1_netKaon_rec_ptSTAR);
    Profile_ptSTAR_Term1C2RecNetCharge->Fill(lV0M, term1_netCharge_rec_ptSTAR);
    Profile_ptSTAR_Term1C2RecNetProton->Fill(lV0M, term1_netProton_rec_ptSTAR);
    Profile_ptSTAR_Term2C2RecNetKaon->Fill(lV0M, term2_netKaon_rec_ptSTAR);
    Profile_ptSTAR_Term2C2RecNetCharge->Fill(lV0M, term2_netCharge_rec_ptSTAR);
    Profile_ptSTAR_Term2C2RecNetProton->Fill(lV0M, term2_netProton_rec_ptSTAR);

    
     
    //++++++++++++++++++++++++++++++++++++++++++++
    //Corrected from Reconstructed
    //++++++++++++++++++++++++++++++++++++++++++++
    
    /*

    //Kaon
    fCorrectedNoKaonPlus_ptmax2 = no_KaonPlus_perevent_ptmax2_corrected;
    fCorrectedNoKaonMinus_ptmax2 = no_KaonMinus_perevent_ptmax2_corrected;
    fCorrectedNoKaonPlus_ptmax3 = no_KaonPlus_perevent_corrected;
    fCorrectedNoKaonMinus_ptmax3 = no_KaonMinus_perevent_corrected;
    //Pion
    fCorrectedNoPionPlus_ptmax2 = no_PionPlus_perevent_ptmax2_corrected;
    fCorrectedNoPionMinus_ptmax2 = no_PionMinus_perevent_ptmax2_corrected;
    fCorrectedNoPionPlus_ptmax3 = no_PionPlus_perevent_corrected;
    fCorrectedNoPionMinus_ptmax3 = no_PionMinus_perevent_corrected;
    //Proton
    fCorrectedNoProtonPlus_ptmax2 = no_ProtonPlus_perevent_ptmax2_corrected;
    fCorrectedNoProtonMinus_ptmax2 = no_ProtonMinus_perevent_ptmax2_corrected;
    fCorrectedNoProtonPlus_ptmax3 = no_ProtonPlus_perevent_corrected;
    fCorrectedNoProtonMinus_ptmax3 = no_ProtonMinus_perevent_corrected;

    //Eff Square factor

    fEffSqrFactrPionMinus_ptmax2 = noByEffSquare_PionMinus_perevent_ptmax2_corrected;
    fEffSqrFactrPionPlus_ptmax2 = noByEffSquare_PionPlus_perevent_ptmax2_corrected;
    fEffSqrFactrProtonMinus_ptmax2 = noByEffSquare_ProtonMinus_perevent_ptmax2_corrected;
    fEffSqrFactrProtonPlus_ptmax2 = noByEffSquare_ProtonPlus_perevent_ptmax2_corrected;
    fEffSqrFactrKaonMinus_ptmax2 = noByEffSquare_KaonMinus_perevent_ptmax2_corrected;
    fEffSqrFactrKaonPlus_ptmax2 = noByEffSquare_KaonPlus_perevent_ptmax2_corrected;

    fEffSqrFactrPionMinus_ptmax3 = noByEffSquare_PionMinus_perevent_corrected;
    fEffSqrFactrPionPlus_ptmax3 = noByEffSquare_PionPlus_perevent_corrected;
    fEffSqrFactrProtonMinus_ptmax3 = noByEffSquare_ProtonMinus_perevent_corrected;
    fEffSqrFactrProtonPlus_ptmax3 = noByEffSquare_ProtonPlus_perevent_corrected;
    fEffSqrFactrKaonMinus_ptmax3 = noByEffSquare_KaonMinus_perevent_corrected;
    fEffSqrFactrKaonPlus_ptmax3 = noByEffSquare_KaonPlus_perevent_corrected;

    */

    Double_t netKaon_corr_ptmax2 = no_KaonPlus_perevent_ptmax2_corrected - no_KaonMinus_perevent_ptmax2_corrected;
    Double_t netCharge_corr_ptmax2 = no_PionPlus_perevent_ptmax2_corrected - no_PionMinus_perevent_ptmax2_corrected + no_KaonPlus_perevent_ptmax2_corrected - no_KaonMinus_perevent_ptmax2_corrected + no_ProtonPlus_perevent_ptmax2_corrected - no_ProtonMinus_perevent_ptmax2_corrected;
    Double_t netProton_corr_ptmax2 = no_ProtonPlus_perevent_ptmax2_corrected - no_ProtonMinus_perevent_ptmax2_corrected;

    Double_t netKaon_corr_ptSTAR = no_KaonPlus_perevent_corrected - no_KaonMinus_perevent_corrected;
    Double_t netCharge_corr_ptSTAR = no_PionPlus_perevent_corrected - no_PionMinus_perevent_corrected + no_KaonPlus_perevent_corrected - no_KaonMinus_perevent_corrected + no_ProtonPlus_perevent_corrected - no_ProtonMinus_perevent_corrected;
    Double_t netProton_corr_ptSTAR = no_ProtonPlus_perevent_corrected - no_ProtonMinus_perevent_corrected;

    Profile_ptmax2_CorrectedNetKaon->Fill(lV0M,netKaon_corr_ptmax2);
    Profile_ptmax2_CorrectedNetCharge->Fill(lV0M,netCharge_corr_ptmax2);
    Profile_ptmax2_CorrectedNetProton->Fill(lV0M,netProton_corr_ptmax2);
    Profile_ptSTAR_CorrectedNetKaon->Fill(lV0M,netKaon_corr_ptSTAR);
    Profile_ptSTAR_CorrectedNetCharge->Fill(lV0M,netCharge_corr_ptSTAR);
    Profile_ptSTAR_CorrectedNetProton->Fill(lV0M,netProton_corr_ptSTAR);

    Double_t term1, term2, term3, term4;

    //net-proton - net-kaon correlation (BS correlation --> <BS>)
    //----------------------------------------------------------------------------------
    term1 = netProton_corr_ptmax2*netKaon_corr_ptmax2;
    Profile_ptmax2_CorrectedNetKaonNetProton->Fill(lV0M, term1);

    term1 = netProton_corr_ptSTAR*netKaon_corr_ptSTAR;
    Profile_ptSTAR_CorrectedNetKaonNetProton->Fill(lV0M, term1);

    //net-proton - net-pion correlation (BQ correlation --> <BQ>)
    //----------------------------------------------------------------------------------
    term1 = netProton_corr_ptmax2*netCharge_corr_ptmax2;
    term2 = no_ProtonPlus_perevent_ptmax2_corrected + no_ProtonMinus_perevent_ptmax2_corrected;
    term3 = noByEffSquare_ProtonPlus_perevent_ptmax2_corrected + noByEffSquare_ProtonMinus_perevent_ptmax2_corrected;
    Profile_ptmax2_CorrectedNetChargeNetProton_term1->Fill(lV0M, term1);
    Profile_ptmax2_CorrectedNetChargeNetProton_term2->Fill(lV0M, term2);
    Profile_ptmax2_CorrectedNetChargeNetProton_term3->Fill(lV0M, term3);

    term1 = netProton_corr_ptSTAR*netCharge_corr_ptSTAR;
    term2 = no_ProtonPlus_perevent_corrected + no_ProtonMinus_perevent_corrected;
    term3 = noByEffSquare_ProtonPlus_perevent_corrected + noByEffSquare_ProtonMinus_perevent_corrected;
    Profile_ptSTAR_CorrectedNetChargeNetProton_term1->Fill(lV0M, term1);
    Profile_ptSTAR_CorrectedNetChargeNetProton_term2->Fill(lV0M, term2);
    Profile_ptSTAR_CorrectedNetChargeNetProton_term3->Fill(lV0M, term3);

    //net-kaon - net-charge correlation (QS correlation --> <QS>)
    //----------------------------------------------------------------------------------
    term1 = netKaon_corr_ptmax2*netCharge_corr_ptmax2;
    term2 = no_KaonPlus_perevent_ptmax2_corrected + no_KaonMinus_perevent_ptmax2_corrected;
    term3 = noByEffSquare_KaonPlus_perevent_ptmax2_corrected + noByEffSquare_KaonMinus_perevent_ptmax2_corrected;
    Profile_ptmax2_CorrectedNetChargeNetKaon_term1->Fill(lV0M, term1);
    Profile_ptmax2_CorrectedNetChargeNetKaon_term2->Fill(lV0M, term2);
    Profile_ptmax2_CorrectedNetChargeNetKaon_term3->Fill(lV0M, term3);

    term1 = netKaon_corr_ptSTAR*netCharge_corr_ptSTAR;
    term2 = no_KaonPlus_perevent_corrected + no_KaonMinus_perevent_corrected;
    term3 = noByEffSquare_KaonPlus_perevent_corrected + noByEffSquare_KaonMinus_perevent_corrected;
    Profile_ptSTAR_CorrectedNetChargeNetKaon_term1->Fill(lV0M, term1);
    Profile_ptSTAR_CorrectedNetChargeNetKaon_term2->Fill(lV0M, term2);
    Profile_ptSTAR_CorrectedNetChargeNetKaon_term3->Fill(lV0M, term3);

    //C2 charge
    //----------------------------------------------------------------------------------
    term1 = TMath::Power(netCharge_corr_ptmax2, 2.0);
    term2 = netCharge_corr_ptmax2;
    term3 = no_PionPlus_perevent_ptmax2_corrected + no_PionMinus_perevent_ptmax2_corrected + no_KaonPlus_perevent_ptmax2_corrected + no_KaonMinus_perevent_ptmax2_corrected + no_ProtonPlus_perevent_ptmax2_corrected + no_ProtonMinus_perevent_ptmax2_corrected;
    term4 = noByEffSquare_ProtonPlus_perevent_ptmax2_corrected + noByEffSquare_ProtonMinus_perevent_ptmax2_corrected + noByEffSquare_PionPlus_perevent_ptmax2_corrected + noByEffSquare_PionMinus_perevent_ptmax2_corrected + noByEffSquare_KaonPlus_perevent_ptmax2_corrected + noByEffSquare_KaonMinus_perevent_ptmax2_corrected;
    Profile_ptmax2_Term1C2CorrectedNetCharge->Fill(lV0M, term1);
    Profile_ptmax2_Term2C2CorrectedNetCharge->Fill(lV0M, term2);
    Profile_ptmax2_Term3C2CorrectedNetCharge->Fill(lV0M, term3);
    Profile_ptmax2_Term4C2CorrectedNetCharge->Fill(lV0M, term4);

    term1 = TMath::Power(netCharge_corr_ptSTAR, 2.0);
    term2 = netCharge_corr_ptSTAR;
    term3 = no_PionPlus_perevent_corrected + no_PionMinus_perevent_corrected + no_KaonPlus_perevent_corrected + no_KaonMinus_perevent_corrected + no_ProtonPlus_perevent_corrected + no_ProtonMinus_perevent_corrected;
    term4 = noByEffSquare_ProtonPlus_perevent_corrected + noByEffSquare_ProtonMinus_perevent_corrected + noByEffSquare_PionPlus_perevent_corrected + noByEffSquare_PionMinus_perevent_corrected + noByEffSquare_KaonPlus_perevent_corrected + noByEffSquare_KaonMinus_perevent_corrected;
    Profile_ptSTAR_Term1C2CorrectedNetCharge->Fill(lV0M, term1);
    Profile_ptSTAR_Term2C2CorrectedNetCharge->Fill(lV0M, term2);
    Profile_ptSTAR_Term3C2CorrectedNetCharge->Fill(lV0M, term3);
    Profile_ptSTAR_Term4C2CorrectedNetCharge->Fill(lV0M, term4);


    //C2 proton
    //----------------------------------------------------------------------------------
    term1 = TMath::Power(netProton_corr_ptmax2, 2.0);
    term2 = netProton_corr_ptmax2;
    term3 = no_ProtonPlus_perevent_ptmax2_corrected + no_ProtonMinus_perevent_ptmax2_corrected;
    term4 = noByEffSquare_ProtonPlus_perevent_ptmax2_corrected + noByEffSquare_ProtonMinus_perevent_ptmax2_corrected;
    Profile_ptmax2_Term1C2CorrectedNetProton->Fill(lV0M, term1);
    Profile_ptmax2_Term2C2CorrectedNetProton->Fill(lV0M, term2);
    Profile_ptmax2_Term3C2CorrectedNetProton->Fill(lV0M, term3);
    Profile_ptmax2_Term4C2CorrectedNetProton->Fill(lV0M, term4);

    term1 = TMath::Power(netProton_corr_ptSTAR, 2.0);
    term2 = netProton_corr_ptSTAR;
    term3 = no_ProtonPlus_perevent_corrected + no_ProtonMinus_perevent_corrected;
    term4 = noByEffSquare_ProtonPlus_perevent_corrected + noByEffSquare_ProtonMinus_perevent_corrected;
    Profile_ptSTAR_Term1C2CorrectedNetProton->Fill(lV0M, term1);
    Profile_ptSTAR_Term2C2CorrectedNetProton->Fill(lV0M, term2);
    Profile_ptSTAR_Term3C2CorrectedNetProton->Fill(lV0M, term3);
    Profile_ptSTAR_Term4C2CorrectedNetProton->Fill(lV0M, term4);

    //C2 kaon
    //----------------------------------------------------------------------------------
    term1 = TMath::Power(netKaon_corr_ptmax2, 2.0);
    term2 = netKaon_corr_ptmax2;
    term3 = no_KaonPlus_perevent_ptmax2_corrected + no_KaonMinus_perevent_ptmax2_corrected;
    term4 = noByEffSquare_KaonPlus_perevent_ptmax2_corrected + noByEffSquare_KaonMinus_perevent_ptmax2_corrected;
    Profile_ptmax2_Term1C2CorrectedNetKaon->Fill(lV0M, term1);
    Profile_ptmax2_Term2C2CorrectedNetKaon->Fill(lV0M, term2);
    Profile_ptmax2_Term3C2CorrectedNetKaon->Fill(lV0M, term3);
    Profile_ptmax2_Term4C2CorrectedNetKaon->Fill(lV0M, term4);

    term1 = TMath::Power(netKaon_corr_ptSTAR, 2.0);
    term2 = netKaon_corr_ptSTAR;
    term3 = no_KaonPlus_perevent_corrected + no_KaonMinus_perevent_corrected;
    term4 = noByEffSquare_KaonPlus_perevent_corrected + noByEffSquare_KaonMinus_perevent_corrected;
    Profile_ptSTAR_Term1C2CorrectedNetKaon->Fill(lV0M, term1);
    Profile_ptSTAR_Term2C2CorrectedNetKaon->Fill(lV0M, term2);
    Profile_ptSTAR_Term3C2CorrectedNetKaon->Fill(lV0M, term3);
    Profile_ptSTAR_Term4C2CorrectedNetKaon->Fill(lV0M, term4);


    //fTreeEvent->Fill();

 
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeList);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPb_v2::GetEvent ()  //event cuts copied from my code written earlier 

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
      return kFALSE;
    }
  
  hNumberOfEvents -> Fill(0.5);

    
  //Standard Event Cuts
  if (!fAODeventCuts.AcceptEvent(fInputEvent)) {
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeList);
    return kFALSE;
  }
  hNumberOfEvents -> Fill(1.5);
  
  
  //Reject Events with Incomplete DAQ
  if (fAODevent->IsIncompleteDAQ())
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeList);
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
      return kFALSE;
    }
  hNumberOfEvents -> Fill(5.5);
  
        
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
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
	return kFALSE;
      }
    hNumberOfEvents -> Fill(7.5);
        
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(8.5);

    
    //SPD Pile-up in Mult Bins
    if (fAODevent->IsPileupFromSPDInMultBins())
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
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
	return kFALSE;
      }
    hNumberOfEvents -> Fill(10.5);
    
        
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(11.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -1.0*fVertexZMax || vertex_tracks->GetZ() > +1.0*fVertexZMax)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
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
	return kFALSE;
      } 
    hNumberOfEvents -> Fill(13.5);
                
    //Selection of Multiplicity Range
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    if (mult_percentile < 0.0 || mult_percentile > 90.0)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeList);
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
Bool_t AliAnalysisTaskCorrPbPb_v2::PassedTrackQualityCuts (AliAODTrack *track)  {
    
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
Bool_t AliAnalysisTaskCorrPbPb_v2::KaonSelector(AliVTrack *track, Double_t nSigmaCut)  {
 
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
Bool_t AliAnalysisTaskCorrPbPb_v2::PionSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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
Bool_t AliAnalysisTaskCorrPbPb_v2::ProtonSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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
Bool_t AliAnalysisTaskCorrPbPb_v2::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
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
Bool_t AliAnalysisTaskCorrPbPb_v2::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kTRUE);
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  return passedTrackPileupCut;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPb_v2::GetMCEffCorrectionHist()
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
 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPb_v2::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
