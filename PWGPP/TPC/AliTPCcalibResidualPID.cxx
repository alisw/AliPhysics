/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>             *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// The task:
// stores TPC PID quantities in a THnSparse
// output can then be used for e.g. dEdx calibration
//
// Author:
// Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//


#include "AliTPCcalibResidualPID.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
//#include "AliStack.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"
#include "AliCentrality.h"
#include "AliAnalysisUtils.h"
#include "THnSparse.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TFile.h"
#include "TSpline.h"
#include "TStyle.h"
#include "AliTPCdEdxInfo.h"
#include "TString.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TVirtualFitter.h"
#include "TTree.h"

using namespace std;

ClassImp(AliTPCcalibResidualPID)

Double_t AliTPCcalibResidualPID::fgCutGeo = 1.;   
Double_t AliTPCcalibResidualPID::fgCutNcr = 0.85; 
Double_t AliTPCcalibResidualPID::fgCutNcl = 0.7;  

//________________________________________________________________________
AliTPCcalibResidualPID::AliTPCcalibResidualPID()
  : AliAnalysisTaskSE(), fESD(0), fMC(0), fOutputContainer(0), fESDtrackCuts(0), fESDtrackCutsV0(0), fPIDResponse(0),
    fNumEtaCorrReqErrorsIssued(0),
    fNumMultCorrReqErrorsIssued(0),
    fUseTPCCutMIGeo(kFALSE),
    fUseMCinfo(kTRUE),
    fIsPbpOrpPb(kFALSE),
    fIsPbPb(kFALSE),
    fZvtxCutEvent(9999.0),
    fV0KineCuts(0x0),
    fAnaUtils(0x0),
    fCutOnProdRadiusForV0el(kTRUE),
    fNumTagsStored(0),
    fV0tags(0x0),
    fV0motherIndex(0x0),
    fV0motherPDG(0x0),
    fProduceAllPadTypes(0),
    fProduceGlobal(0),
    fProduceShortPads(0),
    fProduceMediumPads(0),
    fProduceLongPads(0),
    fProduceOroc(0),
    fHistPidQA(0), 
    fHistPidQAshort(0),
    fHistPidQAmedium(0),
    fHistPidQAlong(0),
    fHistPidQAoroc(0),
    fProduceTPCSignalSparse(0),
    fCorrectdEdxEtaDependence(0),
    fCorrectdEdxMultiplicityDependence(0),
    fCorrectdEdxPileupDependence(kTRUE),
    fThnspTpc(0),
    fWriteAdditionalOutput(kFALSE),
    fQAList(0x0),
    fhInvMassGamma(0x0),
    fhInvMassK0s(0x0),
    fhInvMassLambda(0x0),
    fhInvMassAntiLambda(0x0),
    fhArmenterosAll(0x0),
    fhArmenterosGamma(0x0),
    fhArmenterosK0s(0x0),
    fhArmenterosLambda(0x0),
    fhArmenterosAntiLambda(0x0),
    fHistSharedClusQAV0Pi(0x0),
    fHistSharedClusQAV0Pr(0x0),
    fHistSharedClusQAV0El(0x0),
    fTreeV0El(0x0),
    fTreeV0Pi(0x0),
    fTreeV0Pr(0x0),
    fTree_dEdx_tr(0.),
    fTree_dEdx_nb(0.),
    fTree_dEdx_vs(0.),
    fTree_dEdxExpected_tr(0.),
    fTree_p_TPC_tr(0.),
    fTree_p_TPC_nb(0.),
    fTree_p_TPC_vs(0.),
    fTree_BtimesChargeOverPt_tr(0.),
    fTree_BtimesChargeOverPt_nb(0.),
    fTree_BtimesChargeOverPt_vs(0.),
    fTree_tanTheta_tr(0.),
    fTree_tanTheta_nb(0.),
    fTree_tanTheta_vs(0.),
    fTree_distance_nb(0.),
    fTree_distance_vs(0.)
{
  // default Constructor
   /* fast compilation test
     gSystem->Load("libANALYSIS");
     gSystem->Load("libANALYSISalice");
     .L /lustre/alice/akalweit/train/trunk/util/statsQA/AliTPCcalibResidualPID.cxx++
   */

}


//________________________________________________________________________
AliTPCcalibResidualPID::AliTPCcalibResidualPID(const char *name)
  : AliAnalysisTaskSE(name), fESD(0), fMC(0), fOutputContainer(0), fESDtrackCuts(0), fESDtrackCutsV0(0), fPIDResponse(0),
    fNumEtaCorrReqErrorsIssued(0),
    fNumMultCorrReqErrorsIssued(0),
    fUseTPCCutMIGeo(kFALSE),
    fUseMCinfo(kTRUE),
    fIsPbpOrpPb(kFALSE),
    fIsPbPb(kFALSE),
    fZvtxCutEvent(9999.0),
    fV0KineCuts(0x0),
    fAnaUtils(0x0),
    fCutOnProdRadiusForV0el(kTRUE),
    fNumTagsStored(0),
    fV0tags(0x0),
    fV0motherIndex(0x0),
    fV0motherPDG(0x0),
    fProduceAllPadTypes(0),
    fProduceGlobal(0),
    fProduceShortPads(0),
    fProduceMediumPads(0),
    fProduceLongPads(0),
    fProduceOroc(0),
    fHistPidQA(0),
    fHistPidQAshort(0),
    fHistPidQAmedium(0),
    fHistPidQAlong(0), 
    fHistPidQAoroc(0),
    fProduceTPCSignalSparse(0),
    fCorrectdEdxEtaDependence(0),
    fCorrectdEdxMultiplicityDependence(0),
    fCorrectdEdxPileupDependence(kTRUE),
    fThnspTpc(0),
    fWriteAdditionalOutput(kFALSE),
    fQAList(0x0),
    fhInvMassGamma(0x0),
    fhInvMassK0s(0x0),
    fhInvMassLambda(0x0),
    fhInvMassAntiLambda(0x0),
    fhArmenterosAll(0x0),
    fhArmenterosGamma(0x0),
    fhArmenterosK0s(0x0),
    fhArmenterosLambda(0x0),
    fhArmenterosAntiLambda(0x0),
    fHistSharedClusQAV0Pi(0x0),
    fHistSharedClusQAV0Pr(0x0),
    fHistSharedClusQAV0El(0x0),
    fTreeV0El(0x0),
    fTreeV0Pi(0x0),
    fTreeV0Pr(0x0),
    fTree_dEdx_tr(0.),
    fTree_dEdx_nb(0.),
    fTree_dEdx_vs(0.),
    fTree_dEdxExpected_tr(0.),
    fTree_p_TPC_tr(0.),
    fTree_p_TPC_nb(0.),
    fTree_p_TPC_vs(0.),
    fTree_BtimesChargeOverPt_tr(0.),
    fTree_BtimesChargeOverPt_nb(0.),
    fTree_BtimesChargeOverPt_vs(0.),
    fTree_tanTheta_tr(0.),
    fTree_tanTheta_nb(0.),
    fTree_tanTheta_vs(0.),
    fTree_distance_nb(0.),
    fTree_distance_vs(0.)
{

  //fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
  //fESDtrackCutsV0 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);

  // Constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TObjArray::Class());
  DefineOutput(2, TObjArray::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());

}


//_________________________________________________
AliTPCcalibResidualPID::~AliTPCcalibResidualPID()
{
  delete fOutputContainer;
  fOutputContainer = 0;

  delete fQAList;
  fQAList = 0;

  /*
  delete fTreeV0El;
  fTreeV0El = 0;

  delete fTreeV0Pi;
  fTreeV0Pi = 0;

  delete fTreeV0Pr;
  fTreeV0Pr = 0;
  */

  delete fESDtrackCuts;
  fESDtrackCuts = 0;

  delete fESDtrackCutsV0;
  fESDtrackCutsV0 = 0;

  delete fV0KineCuts;
  fV0KineCuts = 0;
  
  delete fAnaUtils;
  fAnaUtils = 0;

  delete [] fV0tags;
  fV0tags = 0;
  fNumTagsStored = 0;

  delete [] fV0motherIndex;
  fV0motherIndex = 0;

  delete [] fV0motherPDG;
  fV0motherPDG = 0;
}


//________________________________________________________________________
void AliTPCcalibResidualPID::UserCreateOutputObjects()
{
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandler)
    printf("Inputhandler not available \n");
  else
    fPIDResponse = inputHandler->GetPIDResponse();

  // THnSparse binning
  const Int_t kNdim = 9;
  //v0id,  dEdx,  TPCsigele,  TPCsigpion,  TOFbit,  eta, TPCclus, centr, p
  Int_t bins[kNdim] =    {    4,    250,    200,    200,          8,    20,        50,   20,   100};
  Double_t xmin[kNdim] = {  0.5,     30,   -10.,   -10.,       -1.5,   -1.,       60.,   0.,   0.1};
  Double_t xmax[kNdim] = {  4.5,   500.,    10.,    10.,        6.5,    1.,      160.,  100,   4};
  fThnspTpc= new THnSparseF("tpcsignals", "TPC signal;v0id;tpc signal;tpcsigele;tpcsigpion;tofbit;eta;tpcclus;centr;p (GeV/c)", kNdim, bins, xmin, xmax);
  SetAxisNamesFromTitle(fThnspTpc);
  BinLogAxis(fThnspTpc, 8);

  fHistPidQA = InitialisePIDQAHist("fHistPidQA","PID QA", GetIsPbPb());

  fHistPidQAshort  = InitialisePIDQAHist("fHistPidQAshort" ,"PID QA -- short pads", GetIsPbPb());
  fHistPidQAmedium = InitialisePIDQAHist("fHistPidQAmedium","PID QA -- med pads", GetIsPbPb());
  fHistPidQAlong   = InitialisePIDQAHist("fHistPidQAlong"  ,"PID QA -- long pads", GetIsPbPb());
  fHistPidQAoroc   = InitialisePIDQAHist("fHistPidQAoroc"  ,"PID QA -- oroc", GetIsPbPb());

  fOutputContainer = new TObjArray(2);
  fOutputContainer->SetName(GetName());
  fOutputContainer->SetOwner();

  if(fProduceTPCSignalSparse)fOutputContainer->Add(fThnspTpc);
  if(fProduceGlobal)fOutputContainer->Add(fHistPidQA);
  if(fProduceAllPadTypes || fProduceShortPads)fOutputContainer->Add(fHistPidQAshort);
  if(fProduceAllPadTypes || fProduceMediumPads)fOutputContainer->Add(fHistPidQAmedium);
  if(fProduceAllPadTypes || fProduceLongPads)fOutputContainer->Add(fHistPidQAlong);
  if(fProduceAllPadTypes || fProduceOroc)fOutputContainer->Add(fHistPidQAoroc);

  // V0 Kine cuts 
  fV0KineCuts = new AliESDv0KineCuts;
  fV0KineCuts->SetGammaCutChi2NDF(5.);  
  
  fAnaUtils = new AliAnalysisUtils();

  if (fCutOnProdRadiusForV0el) {
    // Only accept V0el with prod. radius within 45 cm -> PID will by systematically biased for larger values!
    Float_t gammaProdVertexRadiusCuts[2] = { 3.0, 45. }; 
    fV0KineCuts->SetGammaCutVertexR(&gammaProdVertexRadiusCuts[0]);
  }


  if (fWriteAdditionalOutput) {
    fQAList = new TObjArray(4);
    fQAList->SetName(GetName());
    fQAList->SetOwner();

    fhInvMassGamma      = new TH1F("fhInvMassGamma", "Invariant Mass of gammas; m_{ee} (GeV/#it{c}^{2}); Entries", 200, 0., 0.2);
    fhInvMassK0s        = new TH1F("fhInvMassK0s", "Invariant Mass of K0s; m_{#pi#pi} (GeV/#it{c}^{2}); Entries;", 200, 0.45, 0.55);
    fhInvMassLambda     = new TH1F("fhInvMassLambda", "Invariant Mass of lambdas; m_{p#pi^{-}} (GeV/#it{c}^{2}); Entries", 200, 1.05, 1.15);
    fhInvMassAntiLambda = new TH1F("fhInvMassAntiLambda", "Invariant Mass of anti-lambdas; m_{#pi^{+}#bar{p}} (GeV/#it{c}^{2}); Entries",
                                  200, 1.05, 1.15);

    fQAList->Add(fhInvMassGamma);
    fQAList->Add(fhInvMassK0s);
    fQAList->Add(fhInvMassLambda);
    fQAList->Add(fhInvMassAntiLambda);


    fhArmenterosAll = new TH2F("fhArmenterosAll",
                              "Armenteros plot all V0s;#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L});#it{q}_{T} (GeV/#it{c})",
                              200, -1., 1., 200, 0., 0.4);
    fhArmenterosGamma = new TH2F("fhArmenterosGamma",
                                "Armenteros plot Gamma;#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L});#it{q}_{T} (GeV/#it{c})",
                                200, -1., 1., 200, 0., 0.4);
    fhArmenterosK0s = new TH2F("fhArmenterosK0s",
                              "Armenteros plot K0s;#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L});#it{q}_{T} (GeV/#it{c})",
                              200, -1., 1., 200, 0., 0.4);
    fhArmenterosLambda = new TH2F("fhArmenterosLambda",
                                "Armenteros plot lambda;#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L});#it{q}_{T} (GeV/#it{c})",
                                200, -1., 1., 200, 0., 0.4);
    fhArmenterosAntiLambda = new TH2F("fhArmenterosAntiLambda",
                                      "Armenteros plot anti-lambda;#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L});#it{q}_{T} (GeV/#it{c})",
                                      200, -1., 1., 200, 0., 0.4);

    fQAList->Add(fhArmenterosAll);
    fQAList->Add(fhArmenterosGamma);
    fQAList->Add(fhArmenterosK0s);
    fQAList->Add(fhArmenterosLambda);
    fQAList->Add(fhArmenterosAntiLambda);


    const Int_t dimQASharedClusters = 4;
    const TString axisTitles[dimQASharedClusters] = { "#it{p}_{TPC} (GeV/#it{c})", "#it{#Delta'}", "#it{N}_{shared cl}", "Pad row"};
    Int_t    binsHistQASharedClusters[dimQASharedClusters] = { 100,  100, 160, 160};
    Double_t xminHistQASharedClusters[dimQASharedClusters] = { 0.3, 0.5,  0,   -1};
    Double_t xmaxHistQASharedClusters[dimQASharedClusters] = { 20,  1.5, 160, 159};

    fHistSharedClusQAV0Pi = new THnSparseF("fHistSharedClusQAV0Pi" ,"PID QA shared clusters - V0 pi", dimQASharedClusters,
                                          binsHistQASharedClusters, xminHistQASharedClusters, xmaxHistQASharedClusters);
    BinLogAxis(fHistSharedClusQAV0Pi, 0);
    for (Int_t i = 0; i < dimQASharedClusters; i++)
      fHistSharedClusQAV0Pi->GetAxis(i)->SetTitle(axisTitles[i].Data());
    fQAList->Add(fHistSharedClusQAV0Pi);

    fHistSharedClusQAV0Pr = new THnSparseF("fHistSharedClusQAV0Pr" ,"PID QA shared clusters - V0 pr", dimQASharedClusters,
                                          binsHistQASharedClusters, xminHistQASharedClusters, xmaxHistQASharedClusters);
    BinLogAxis(fHistSharedClusQAV0Pr, 0);
    for (Int_t i = 0; i < dimQASharedClusters; i++)
      fHistSharedClusQAV0Pi->GetAxis(i)->SetTitle(axisTitles[i].Data());
    fQAList->Add(fHistSharedClusQAV0Pr);

    fHistSharedClusQAV0El = new THnSparseF("fHistSharedClusQAV0El" ,"PID QA shared clusters - V0 el", dimQASharedClusters,
                                          binsHistQASharedClusters, xminHistQASharedClusters, xmaxHistQASharedClusters);
    BinLogAxis(fHistSharedClusQAV0El, 0);
    for (Int_t i = 0; i < dimQASharedClusters; i++)
      fHistSharedClusQAV0Pi->GetAxis(i)->SetTitle(axisTitles[i].Data());
    fQAList->Add(fHistSharedClusQAV0El);

    OpenFile(3);
    fTreeV0El = new TTree("treeV0El", "V0 el together with info from closest neighbour");
    fTreeV0El->Branch("dEdx_tr", &fTree_dEdx_tr);
    fTreeV0El->Branch("dEdx_nb", &fTree_dEdx_nb);
    fTreeV0El->Branch("dEdx_vs", &fTree_dEdx_vs);
    fTreeV0El->Branch("dEdxExpected_tr", &fTree_dEdxExpected_tr);
    fTreeV0El->Branch("p_TPC_tr", &fTree_p_TPC_tr);
    fTreeV0El->Branch("p_TPC_nb", &fTree_p_TPC_nb);
    fTreeV0El->Branch("p_TPC_vs", &fTree_p_TPC_vs);
    fTreeV0El->Branch("BtimesChargeOverPt_tr", &fTree_BtimesChargeOverPt_tr);
    fTreeV0El->Branch("BtimesChargeOverPt_nb", &fTree_BtimesChargeOverPt_nb);
    fTreeV0El->Branch("BtimesChargeOverPt_vs", &fTree_BtimesChargeOverPt_vs);
    fTreeV0El->Branch("tanTheta_tr", &fTree_tanTheta_tr);
    fTreeV0El->Branch("tanTheta_nb", &fTree_tanTheta_nb);
    fTreeV0El->Branch("tanTheta_vs", &fTree_tanTheta_vs);
    fTreeV0El->Branch("distance_nb", &fTree_distance_nb);
    fTreeV0El->Branch("distance_vs", &fTree_distance_vs);

    OpenFile(4);
    fTreeV0Pi = new TTree("treeV0Pi", "V0 pi together with info from closest neighbour");
    fTreeV0Pi->Branch("dEdx_tr", &fTree_dEdx_tr);
    fTreeV0Pi->Branch("dEdx_nb", &fTree_dEdx_nb);
    fTreeV0Pi->Branch("dEdx_vs", &fTree_dEdx_vs);
    fTreeV0Pi->Branch("dEdxExpected_tr", &fTree_dEdxExpected_tr);
    fTreeV0Pi->Branch("p_TPC_tr", &fTree_p_TPC_tr);
    fTreeV0Pi->Branch("p_TPC_nb", &fTree_p_TPC_nb);
    fTreeV0Pi->Branch("p_TPC_vs", &fTree_p_TPC_vs);
    fTreeV0Pi->Branch("BtimesChargeOverPt_tr", &fTree_BtimesChargeOverPt_tr);
    fTreeV0Pi->Branch("BtimesChargeOverPt_nb", &fTree_BtimesChargeOverPt_nb);
    fTreeV0Pi->Branch("BtimesChargeOverPt_vs", &fTree_BtimesChargeOverPt_vs);
    fTreeV0Pi->Branch("tanTheta_tr", &fTree_tanTheta_tr);
    fTreeV0Pi->Branch("tanTheta_nb", &fTree_tanTheta_nb);
    fTreeV0Pi->Branch("tanTheta_vs", &fTree_tanTheta_vs);
    fTreeV0Pi->Branch("distance_nb", &fTree_distance_nb);
    fTreeV0Pi->Branch("distance_vs", &fTree_distance_vs);

    OpenFile(5);
    fTreeV0Pr = new TTree("treeV0Pr", "V0 pr together with info from closest neighbour");
    fTreeV0Pr->Branch("dEdx_tr", &fTree_dEdx_tr);
    fTreeV0Pr->Branch("dEdx_nb", &fTree_dEdx_nb);
    fTreeV0Pr->Branch("dEdx_vs", &fTree_dEdx_vs);
    fTreeV0Pr->Branch("dEdxExpected_tr", &fTree_dEdxExpected_tr);
    fTreeV0Pr->Branch("p_TPC_tr", &fTree_p_TPC_tr);
    fTreeV0Pr->Branch("p_TPC_nb", &fTree_p_TPC_nb);
    fTreeV0Pr->Branch("p_TPC_vs", &fTree_p_TPC_vs);
    fTreeV0Pr->Branch("BtimesChargeOverPt_tr", &fTree_BtimesChargeOverPt_tr);
    fTreeV0Pr->Branch("BtimesChargeOverPt_nb", &fTree_BtimesChargeOverPt_nb);
    fTreeV0Pr->Branch("BtimesChargeOverPt_vs", &fTree_BtimesChargeOverPt_vs);
    fTreeV0Pr->Branch("tanTheta_tr", &fTree_tanTheta_tr);
    fTreeV0Pr->Branch("tanTheta_nb", &fTree_tanTheta_nb);
    fTreeV0Pr->Branch("tanTheta_vs", &fTree_tanTheta_vs);
    fTreeV0Pr->Branch("distance_nb", &fTree_distance_nb);
    fTreeV0Pr->Branch("distance_vs", &fTree_distance_vs);
  }

  PostData(1,fOutputContainer);
  if (fWriteAdditionalOutput) {
    PostData(2,fQAList);
    PostData(3,fTreeV0El);
    PostData(4,fTreeV0Pi);
    PostData(5,fTreeV0Pr);
  }
}


//_____________________________________________________________________________
void AliTPCcalibResidualPID::UserExec(Option_t *)
{
    //calls the Process function
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      printf("ERROR: Could not get ESDInputHandler \n");
    }
    else fESD = (AliESDEvent*)esdH->GetEvent();

    // If MC, forward MCevent
    fMC = dynamic_cast<AliMCEvent*>(MCEvent());
    //
    Process(fESD, fMC);
    //
    PostData(1,fOutputContainer);
    if (fWriteAdditionalOutput) {
      PostData(2,fQAList);
      PostData(3,fTreeV0El);
      PostData(4,fTreeV0Pi);
      PostData(5,fTreeV0Pr);
    }
}


//________________________________________________________________________
THnSparseF* AliTPCcalibResidualPID::InitialisePIDQAHist(TString name, TString title, Bool_t IsPbPb)
{
  // Initialise a pidQA histo

  //
  // 0.ptot, 1.tpcSig, 2.particle ID, 3. assumed particle, 4. nSigmaTPC (4x), 5. nSigmaTOF (4x), 6. centrality
  // Concerning 2 (particle ID):
  // - in case of MC: Contains MC ID. Bin 1-4 (<=> Slot 0-3): el, pi, ka, pr
  // - in case of data: Contains V0 particles + bin with all others: Bin 1-4 (<=> Slot 0-3): All non-V0s, V0-el, V0-pi, V0-pr
  //
  title.Append(";p (GeV/c);tpc signal;particle ID;assumed particle;nSigmaTPC;nSigmaTOF;centrality");
  const Int_t kNdim = 7;
  Int_t    binsHistQA[kNdim] = {135, 1980,    4,    5, 40, 10,  IsPbPb ? 40 : 40 };
  Double_t xminHistQA[kNdim] = {0.1,   20., -0.5, -0.5, -10., -5.,   0.};
  Double_t xmaxHistQA[kNdim] = {50., 2000.,  3.5,  4.5,  10.,  5., IsPbPb ? 20000. : 4000.};
  THnSparseF* h = new THnSparseF(name.Data(), title.Data(), kNdim, binsHistQA, xminHistQA, xmaxHistQA);
  BinLogAxis(h, 0);
  SetAxisNamesFromTitle(h);

  return h;
}


//________________________________________________________________________
void AliTPCcalibResidualPID::Process(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
{

  //called for each event
  if (!esdEvent) {
    Printf("ERROR: esdEvent not available"); 
    return;
  }

  if (!fPIDResponse || !fV0KineCuts || !fAnaUtils) {
    Printf("ERROR: No PIDresponse, v0KineCuts or AliAnalysisUtils!");
    return;
  }
  
  if (fAnaUtils->IsPileUpSPD(esdEvent))
    return;

  Float_t centralityFper=99;

  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityFper = esdCentrality->GetCentralityPercentile("V0M");

  if (!GetVertexIsOk(esdEvent))
    return;

  const AliESDVertex* fESDvertex = esdEvent->GetPrimaryVertexTracks(); 
  if (!fESDvertex)
    return;

  Int_t ncontr = fESDvertex->GetNContributors();

  if (ncontr <= 0)
    return;

  // Fill V0 arrays with V0s
  FillV0PIDlist(esdEvent);

  // Array with flags wheter QA for this V0 was already done or not
  const Int_t numV0s = esdEvent->GetNumberOfV0s();
  Bool_t v0QAadded[numV0s];
  for (Int_t i = 0; i < numV0s; i++)
    v0QAadded[i] = kFALSE;

  Int_t nTotTracks = esdEvent->GetNumberOfTracks();
  const Int_t nTotESDTracks = esdEvent->GetNumberOfESDTracks();

  const Double_t magField = esdEvent->GetMagneticField();

  const Bool_t etaCorrAvail    = fPIDResponse->UseTPCEtaCorrection();
  const Bool_t multCorrAvail   = fPIDResponse->UseTPCMultiplicityCorrection();
  const Bool_t pileupCorrAvail = fPIDResponse->UseTPCPileupCorrection();

  for (Int_t iTracks = 0; iTracks < nTotTracks; iTracks++){//begin track loop 
    AliESDtrack *trackESD = esdEvent->GetTrack(iTracks);
    if(!trackESD) {
      Printf("ERROR: Could not receive track %d (esd loop)", iTracks);
      continue;
    }
    if((TMath::Abs(trackESD->Eta())) > 0.9)
      continue;

     // Do not distinguish positively and negatively charged V0's
    Char_t v0tagAllCharges = TMath::Abs(GetV0tag(iTracks));
    if (v0tagAllCharges == -99) {
      AliError(Form("Problem with V0 tag list (requested V0 track for track %d from %d (list status %d))!", iTracks, esdEvent->GetNumberOfTracks(),
                    fV0tags != 0x0));
      continue;
    }

    Bool_t isV0el = v0tagAllCharges == 14;
    Bool_t isV0pi = v0tagAllCharges == 15;
    Bool_t isV0pr = v0tagAllCharges == 16;
    Bool_t isV0 = isV0el || isV0pi || isV0pr;

    if (mcEvent && fUseMCinfo) {
      // For MC, do not look for V0s, i.e. demand the non-V0 track cuts
      if (fESDtrackCuts && !fESDtrackCuts->AcceptTrack(trackESD)) continue;

      if (fUseTPCCutMIGeo) {
        // If cut on geometry is active, don't cut on number of clusters, since such a cut is already included.
        if (!TPCCutMIGeo(trackESD, esdEvent))
          continue;
      }
      else {
        // If cut on geometry is not active, always cut on num clusters
        if (trackESD->GetTPCsignalN() < 60)
          continue;
      }
    }
    else {
      // For data, take V0 AND non-V0 with separate cuts
      //if (!isV0 && fESDtrackCuts && !fESDtrackCuts->AcceptTrack(trackESD)) continue;
      if (isV0) {
        if (fESDtrackCutsV0 && !fESDtrackCutsV0->AcceptTrack(trackESD))
          continue;
      }
      else {
        if (fESDtrackCuts && !fESDtrackCuts->AcceptTrack(trackESD))
          continue;
      }

      if (fUseTPCCutMIGeo) {
        // If cut on geometry is active, don't cut on number of clusters, since such a cut is already included.
        if (!isV0) {
          if (!TPCCutMIGeo(trackESD, esdEvent))
            continue;
        }
        else {
          // One should not cut on geometry for V0's since they have different topology. (Loosely) cut on num of clusters instead.
          if (trackESD->GetTPCsignalN() < 60)
            continue;
        }
      }
      else {
        // If cut on geometry is not active, always cut on num clusters
        if (trackESD->GetTPCsignalN() < 60)
          continue;
      }
    }

    const AliExternalTrackParam *paramIn = trackESD->GetInnerParam();
    Float_t precin=-1;
    if(paramIn) {
      precin=paramIn->GetP();
    } else {
      printf("Inner param not found for track %d - skipping!\n", iTracks);
      continue;
    }
    Int_t precdefault=CompareFloat(precin,-1);
    if(precdefault==1) continue; // momentum cut
    //

    AliMCParticle* mcTrack = 0x0;
    Int_t particleID = -1;
    Int_t v0id = 0;

    if (mcEvent && fUseMCinfo) {
      // MC - particle ID = MC ID
      Int_t label = trackESD->GetLabel();

      if (label < 0)
        continue; // If MC is available, reject tracks with negative ESD label

      mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(TMath::Abs(label)));
      if (!mcTrack) {
        Printf("ERROR: Could not receive mcTrack with label %d for track %d", label, iTracks);
        continue;
      }

      /*// Only accept MC primaries
      if (!mcEvent->Stack()->IsPhysicalPrimary(TMath::Abs(label))) {
        continue;
      }*/
      
      Int_t pdgAbs = TMath::Abs(mcTrack->PdgCode());
      
      if (pdgAbs == 11) { // electron
        particleID = 0;
      }
      else if (pdgAbs == 211) { // pion
        particleID = 1;
      }
      else if (pdgAbs == 321) { // kaon
        particleID = 2;
      }
      else if (pdgAbs == 2212) { // proton
        particleID = 3;
      }
      else
        continue; // Reject all other particles
    }
    else {
      // Data - particle ID = V0 ID (if V0)
      if (isV0pi) { // pion
        particleID = 2;
        v0id = 2;
      }
      else if (isV0el) { // electron
        particleID = 1;
        v0id = 1;
      }
      else if (isV0pr) { // proton
        particleID = 3;
        v0id = 3;
      }
      else { // No V0-ID available (species must be determined via TPC, TOF, ...)
        particleID = 0;
      }
    }


//     Float_t tpcsignal=trackESD->GetTPCsignal();
    Float_t tpcsignal=fPIDResponse->GetTPCResponse().GetTrackdEdx(trackESD);
    Int_t tofbit=-1;
    Int_t iTOFpid=0;
    Int_t ikTIME=0;
    Double_t tpcnsigmaele=-10;
    Double_t tpcnsigmapion=-10;
    //
    if ((trackESD->GetStatus() & AliESDtrack::kTOFpid)) iTOFpid = 1;
    if ((trackESD->GetStatus() & AliESDtrack::kTIME)) ikTIME = 1;
    Float_t time0 = fPIDResponse->GetTOFResponse().GetTimeZero();
    //
    if((iTOFpid==1) &&(ikTIME==1)){
      Double_t tofnsigmaele= fPIDResponse->NumberOfSigmasTOF(trackESD,AliPID::kElectron, time0);
      Double_t tofnsigmapion=fPIDResponse->NumberOfSigmasTOF(trackESD,AliPID::kPion, time0);
      if (TMath::Abs(tofnsigmapion)<3) tofbit = 1;
      if (TMath::Abs(tofnsigmaele)<3)  tofbit = 0;
      tpcnsigmaele=fPIDResponse->NumberOfSigmasTPC(trackESD,AliPID::kElectron);
      tpcnsigmapion=fPIDResponse->NumberOfSigmasTPC(trackESD,AliPID::kPion);
    }
    //
    Int_t tpcnclusN=trackESD->GetTPCsignalN();
    Double_t eta=trackESD->Eta();
    
    //
    if(fProduceTPCSignalSparse){
      Double_t contentSignal[9];
      contentSignal[0]=v0id;
      contentSignal[1]=tpcsignal;
      contentSignal[2]=tpcnsigmaele;
      contentSignal[3]=tpcnsigmapion;
      contentSignal[4]=tofbit;
      contentSignal[5]=eta;
      contentSignal[6]=tpcnclusN;
      contentSignal[7]=centralityFper;
      contentSignal[8]=precin;
      //
      if (fThnspTpc->GetEntries() < 1e6) fThnspTpc->Fill(contentSignal);
    }//should it be created or not



    //
    // 2nd THnSparse
    //
    Double_t tpcQA[5] = {fPIDResponse->NumberOfSigmasTPC(trackESD, AliPID::kElectron),
                         fPIDResponse->NumberOfSigmasTPC(trackESD, AliPID::kPion),
                         fPIDResponse->NumberOfSigmasTPC(trackESD, AliPID::kKaon),
                         fPIDResponse->NumberOfSigmasTPC(trackESD, AliPID::kProton),
                         0};
    Double_t tofQA[5] = {fPIDResponse->NumberOfSigmasTOF(trackESD, AliPID::kElectron, time0),
                         fPIDResponse->NumberOfSigmasTOF(trackESD, AliPID::kPion, time0),
                         fPIDResponse->NumberOfSigmasTOF(trackESD, AliPID::kKaon, time0),
                         fPIDResponse->NumberOfSigmasTOF(trackESD, AliPID::kProton, time0),
                         0 };
       
    // id 5 is just again Kaons in restricted eta range
    tpcQA[4] = tpcQA[2];
    tofQA[4] = tofQA[2];
    
    //
    // dE/dx in different pad regions
    //
    AliTPCdEdxInfo * infoTpcPid = trackESD->GetTPCdEdxInfo();
    Double32_t signal[4]; Char_t ncl[3]; Char_t nrows[3];
    if (infoTpcPid) {
      infoTpcPid->GetTPCSignalRegionInfo(signal, ncl, nrows);
    } else {
      for(Int_t iarr = 0; iarr < 3; iarr++) {
        signal[iarr] = 0;
        ncl[iarr] = 0;
        nrows[iarr] = 0;
      }
      signal[3] = 0;
    }
    //
    UInt_t status = trackESD->GetStatus();
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout && hasTOFtime) hasTOF = kTRUE;
    Float_t length = trackESD->GetIntegratedLength();
    // Length check only for primaries!
    if (length < 350 && !isV0) hasTOF = kFALSE;
    //
    
    if (!hasTOF) {
      // Make sure that number of sigmas is large in this case, so that the track will be rejected if a TOF cut is applied
      for (Int_t i = 0; i < 5; i++) {
        tofQA[i] = 999;
      }
    }


    Double_t processedTPCsignal[5] = { tpcsignal, tpcsignal, tpcsignal, tpcsignal, tpcsignal };
  
    if (fCorrectdEdxEtaDependence && fNumEtaCorrReqErrorsIssued < 23 && !etaCorrAvail) {
      AliError("TPC eta correction requested, but was not activated in PID response (most likely not available)!");
      fNumEtaCorrReqErrorsIssued++;
      if (fNumEtaCorrReqErrorsIssued == 23)
        AliError("Ignoring this error from now on!");
    }
    
    if (fCorrectdEdxMultiplicityDependence && fNumMultCorrReqErrorsIssued < 23 && !multCorrAvail) {
      AliError("TPC multiplicity correction requested, but was not activated in PID response (most likely not available)!");
      fNumMultCorrReqErrorsIssued++;
      if (fNumMultCorrReqErrorsIssued == 23)
        AliError("Ignoring this error from now on!");
    }
    
    AliTPCPIDResponse& tpcResponse = fPIDResponse->GetTPCResponse();
    const Bool_t etaCorrected = fCorrectdEdxEtaDependence && etaCorrAvail;
    const Bool_t multCorrected = fCorrectdEdxMultiplicityDependence && multCorrAvail;
    const Bool_t pileupCorrected = fCorrectdEdxPileupDependence && pileupCorrAvail;

    processedTPCsignal[0] = tpcResponse.GetCorrectedTrackdEdx(trackESD, AliPID::kElectron, etaCorrected, multCorrected, pileupCorrected);
    processedTPCsignal[1] = tpcResponse.GetCorrectedTrackdEdx(trackESD, AliPID::kPion,     etaCorrected, multCorrected, pileupCorrected);
    processedTPCsignal[2] = tpcResponse.GetCorrectedTrackdEdx(trackESD, AliPID::kKaon,     etaCorrected, multCorrected, pileupCorrected);
    processedTPCsignal[3] = tpcResponse.GetCorrectedTrackdEdx(trackESD, AliPID::kProton,   etaCorrected, multCorrected, pileupCorrected);

    // id 5 is just again Kaons in restricted eta range
    processedTPCsignal[4] = processedTPCsignal[2];
    
    for(Int_t iPart = 0; iPart < 5; iPart++) {
      // Only accept "Kaons" within |eta| < 0.2 for index 4 in case of data (no contamination in case of MC, so this index is not used)
      if (iPart == 4 && ((mcEvent && fUseMCinfo) || abs(trackESD->Eta()) > 0.2)) {
        continue;
      }
      
      Double_t vecHistQA[7] = {precin, processedTPCsignal[iPart], (Double_t)particleID, (Double_t)iPart, tpcQA[iPart], tofQA[iPart],
                               (Double_t)nTotESDTracks};
      if (fProduceGlobal) fHistPidQA->Fill(vecHistQA);
      vecHistQA[1] = signal[0]; vecHistQA[4] = ncl[0];
      if ((fProduceAllPadTypes || fProduceShortPads)) fHistPidQAshort->Fill(vecHistQA);
      //
      vecHistQA[1] = signal[1]; vecHistQA[4] = ncl[1];
      if ((fProduceAllPadTypes || fProduceMediumPads)) fHistPidQAmedium->Fill(vecHistQA);
      //
      vecHistQA[1] = signal[2]; vecHistQA[4] = ncl[2];
      if ((fProduceAllPadTypes || fProduceLongPads)) fHistPidQAlong->Fill(vecHistQA);
      //
      vecHistQA[1] = signal[3]; vecHistQA[4] = nrows[1] + nrows[2];
      if ((fProduceAllPadTypes || fProduceOroc)) fHistPidQAoroc->Fill(vecHistQA);      
    }
    
    if (fWriteAdditionalOutput && isV0) {
      
      // Find the closest neighbour and fill the information into the trees.
      fTree_p_TPC_tr = precin;
      fTree_dEdx_tr = tpcsignal;
      fTree_BtimesChargeOverPt_tr = TMath::Abs(magField) * trackESD->GetSigned1Pt();
      fTree_tanTheta_tr = trackESD->GetInnerParam()->GetTgl();
      
      Int_t closestNeighbourIndex = -1;
      const Double_t tpcInnerRadius = 85.;
      const Int_t parIndexGlobalPhi = 8;
      const Double_t z_tr = trackESD->GetInnerParam()->GetZ();
      const Double_t phi_tr = trackESD->GetInnerParam()->GetParameterAtRadius(tpcInnerRadius, magField, parIndexGlobalPhi);
      
      const Double_t distanceCut = 12.; // Cut on distance between track to closest neighbour
      
      Double_t z_nb = 9999;
      Double_t phi_nb = 9999;
      Double_t delta_z = 9999;
      Double_t delta_phi = 9999;
      
      fTree_distance_nb = -9999;
      
      for (Int_t jTracks = 0; jTracks < nTotTracks; jTracks++) {
        if (iTracks == jTracks)
          continue; // Exclude track we are looking at right now from neighbour list
        
        AliESDtrack *track_nb= esdEvent->GetTrack(jTracks);
        if(!track_nb)
          continue;
        
        const AliExternalTrackParam *paramIn_nb = track_nb->GetInnerParam();
        if (!paramIn_nb)
          continue;
        
        z_nb =  paramIn_nb->GetZ();
        delta_z = TMath::Abs(z_nb - z_tr);
        
        if (delta_z >= distanceCut)
          continue;
        
        phi_nb = paramIn_nb->GetParameterAtRadius(tpcInnerRadius, magField, parIndexGlobalPhi);
        delta_phi = TMath::Abs(phi_nb - phi_tr);
        
        const Double_t delta_r_phi = tpcInnerRadius * delta_phi;
        
        if (delta_r_phi >= distanceCut)
          continue;
        
        
        const Double_t distance = TMath::Sqrt(delta_z * delta_z + delta_r_phi * delta_r_phi);
        
        if (distance >= distanceCut)
          continue;
        
        if (closestNeighbourIndex < 0 || distance < fTree_distance_nb) {
          fTree_distance_nb = distance;
          closestNeighbourIndex = jTracks;
        }
      }
      
      if (closestNeighbourIndex >= 0 && closestNeighbourIndex < nTotTracks) {
        AliESDtrack *track_nb= esdEvent->GetTrack(closestNeighbourIndex);
//         fTree_dEdx_nb = track_nb->GetTPCsignal();
        fTree_dEdx_nb = fPIDResponse->GetTPCResponse().GetTrackdEdx(track_nb);
        fTree_p_TPC_nb = track_nb->GetInnerParam()->GetP();
        fTree_BtimesChargeOverPt_nb = TMath::Abs(magField) * track_nb->GetSigned1Pt();
        fTree_tanTheta_nb = track_nb->GetInnerParam()->GetTgl();
      }
      else {
        // No neighbour in this range -> nevertheless store the track, but set neighbour values to invalid
        fTree_dEdx_nb = -999.;
        fTree_p_TPC_nb = -999.;
        fTree_BtimesChargeOverPt_nb = -999.;
        fTree_tanTheta_nb = -999.;
        fTree_distance_nb = -999.;
      }
      
      
      // Find the other V0 daughter and store its values
      const Int_t motherIndex = GetV0motherIndex(iTracks);
      AliESDv0* v0Mother = (AliESDv0*)esdEvent->GetV0(motherIndex);
      if (!v0Mother) {
        printf("Error: Track tagged as V0 daughter, but V0 mother not found!\n");
        continue;
      }
      
      const Int_t iTrackP = v0Mother->GetPindex();  // positive track
      const Int_t iTrackN = v0Mother->GetNindex();  // negative track_nb
      
      Int_t iTrackV0sister = -1;
      if (iTrackP == iTracks)
        iTrackV0sister = iTrackN;
      else if (iTrackN == iTracks)
        iTrackV0sister = iTrackP;
      else {
        printf("Error: V0 sister relations are wrong!\n");
        continue;
      }
      
      fTree_dEdx_vs = -999.;
      fTree_p_TPC_vs = -999.;
      fTree_BtimesChargeOverPt_vs = -999.;
      fTree_tanTheta_vs = -999.;
      fTree_distance_vs = -999.;
        
      AliESDtrack *track_vs = esdEvent->GetTrack(iTrackV0sister);
      if (track_vs) {
        const AliExternalTrackParam *paramIn_vs = track_vs->GetInnerParam();
        if (paramIn_vs) {
//           fTree_dEdx_vs = track_vs->GetTPCsignal();
          fTree_dEdx_vs = fPIDResponse->GetTPCResponse().GetTrackdEdx(track_vs);
          fTree_p_TPC_vs = paramIn_vs->GetP();
          fTree_BtimesChargeOverPt_vs = TMath::Abs(magField) * track_vs->GetSigned1Pt();
          fTree_tanTheta_vs = paramIn_vs->GetTgl();
          
          // Calculate distance
          const Double_t z_vs =  paramIn_vs->GetZ();
          const Double_t delta_z_vs = TMath::Abs(z_vs - z_tr);
          
          const Double_t phi_vs = paramIn_vs->GetParameterAtRadius(tpcInnerRadius, magField, parIndexGlobalPhi);
          const Double_t delta_phi_vs = TMath::Abs(phi_vs - phi_tr);
          
          const Double_t delta_r_phi_vs = tpcInnerRadius * delta_phi_vs;
          
          fTree_distance_vs = TMath::Sqrt(delta_z_vs * delta_z_vs + delta_r_phi_vs * delta_r_phi_vs);
        }
      }
      
      if (isV0el) {
        fTree_dEdxExpected_tr = fPIDResponse->GetTPCResponse().GetExpectedSignal(trackESD, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, 
                                                                                 fPIDResponse->UseTPCEtaCorrection(),
                                                                                 fPIDResponse->UseTPCMultiplicityCorrection());
        if (fTree_dEdxExpected_tr > 0)
          fTreeV0El->Fill();
      }
      else if (isV0pi) {
        fTree_dEdxExpected_tr = fPIDResponse->GetTPCResponse().GetExpectedSignal(trackESD, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, 
                                                                                 fPIDResponse->UseTPCEtaCorrection(),
                                                                                 fPIDResponse->UseTPCMultiplicityCorrection());
        if (fTree_dEdxExpected_tr > 0)
          fTreeV0Pi->Fill();
      }
      else if (isV0pr) {
        fTree_dEdxExpected_tr = fPIDResponse->GetTPCResponse().GetExpectedSignal(trackESD, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, 
                                                                                 fPIDResponse->UseTPCEtaCorrection(),
                                                                                 fPIDResponse->UseTPCMultiplicityCorrection());
        if (fTree_dEdxExpected_tr > 0)
          fTreeV0Pr->Fill();
      }
      
      // Fill QA hists for shared clusters dE/dx
      Double_t vecHistQA[4] = { precin, -1., (Double_t)trackESD->GetTPCSharedMap().CountBits(), -1. };
      THnSparse* currHist = fHistSharedClusQAV0Pi;
      
      if (isV0pi) {
        Double_t expectedDeDx = fPIDResponse->GetTPCResponse().GetExpectedSignal(trackESD, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault,
                                                                                 fPIDResponse->UseTPCEtaCorrection(),
                                                                                 fPIDResponse->UseTPCMultiplicityCorrection()); 
        if (expectedDeDx > 0 ) {
          vecHistQA[1] = tpcsignal / expectedDeDx;
          currHist = fHistSharedClusQAV0Pi;
        }
      }
      else if (isV0pr) {
        Double_t expectedDeDx = fPIDResponse->GetTPCResponse().GetExpectedSignal(trackESD, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,
                                                                                 fPIDResponse->UseTPCEtaCorrection(),
                                                                                 fPIDResponse->UseTPCMultiplicityCorrection()); 
        if (expectedDeDx > 0 ) {
          vecHistQA[1] = tpcsignal / expectedDeDx;
          currHist = fHistSharedClusQAV0Pr;
        }
      }
      else if (isV0el) {
        Double_t expectedDeDx = fPIDResponse->GetTPCResponse().GetExpectedSignal(trackESD, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,
                                                                                 fPIDResponse->UseTPCEtaCorrection(), 
                                                                                 fPIDResponse->UseTPCMultiplicityCorrection()); 
        if (expectedDeDx > 0 ) {
          vecHistQA[1] = tpcsignal / expectedDeDx;
          currHist = fHistSharedClusQAV0El;
        }
      }
      
      Int_t iRowInd = -1;
      // iRowInd == -1 for "all rows w/o multiple counting"
      currHist->Fill(vecHistQA);
      
      // Fill hist for every pad row with shared cluster
      for (iRowInd = 0; iRowInd < 159; iRowInd++) {
        if (trackESD->GetTPCSharedMap().TestBitNumber(iRowInd)) {
          vecHistQA[3] = iRowInd;
          currHist->Fill(vecHistQA);
        }
      }
      
      
      
      // Check, whether the QA for this V0 mother was already filled into the histograms (is the case if the other daughter has already
      // been processed). If not, set flag to kTRUE and fill the QA.
      const Int_t iV0 = GetV0motherIndex(iTracks);
      if (!v0QAadded[iV0]) {
        v0QAadded[iV0] = kTRUE;
      
        AliESDv0* esdV0 = (AliESDv0*)esdEvent->GetV0(iV0);
        
        if (!esdV0) {
          printf("Error: V0 tagged, but does not exist!\n");
        }
        else {
          Float_t armVar[2] = {0.0,0.0};
          fV0KineCuts->Armenteros(esdV0, armVar);
          
          const Int_t motherPDG = GetV0motherPDG(iTracks);
          
          // Mother and daughter match the requirements, otherwise it wouldn't be tagged as a V0. Now just fill the QA histos.
          fhArmenterosAll->Fill(armVar[0], armVar[1]);
          //if (esdV0->TestBit(BIT(14))) {
          if (TMath::Abs(motherPDG) == 22) {
            fhInvMassGamma->Fill(esdV0->GetEffMass(AliPID::kElectron, AliPID::kElectron));
            fhArmenterosGamma->Fill(armVar[0], armVar[1]);
          }
          else if (TMath::Abs(motherPDG) == 310) {
          //else if (esdV0->TestBit(BIT(15))) {
            fhInvMassK0s->Fill(esdV0->GetEffMass(AliPID::kPion, AliPID::kPion));
            fhArmenterosK0s->Fill(armVar[0], armVar[1]);
          }
          else if (motherPDG == 3122) {
          //else if (esdV0->TestBit(BIT(16))) {
            fhInvMassLambda->Fill(esdV0->GetEffMass(AliPID::kProton, AliPID::kPion));
            fhArmenterosLambda->Fill(armVar[0], armVar[1]);
          }
          else if (motherPDG == -3122) {
          //else if (esdV0->TestBit(BIT(17))) {
            fhInvMassAntiLambda->Fill(esdV0->GetEffMass(AliPID::kPion, AliPID::kProton));
            fhArmenterosAntiLambda->Fill(armVar[0], armVar[1]);
          }
        }
      }
    }
  } //end track loop 
  

  // Clear the V0 PID arrays
  ClearV0PIDlist();
}      


//________________________________________________________________________
Bool_t AliTPCcalibResidualPID::ProcessV0Tree(TTree* tree, THnSparseF* h, const Int_t recoPass/*=4*/,
                                             const TString runList/*=""*/, const Bool_t excludeRuns/*=kFALSE*/)
{
  // Fill tracks from tree into h. Returns kTRUE on success.
  // recoPass is required to propbely initialise the AliPIDResponse
  // runList can either be a whith or black list of runs, ignored if empty
  // if excludeRuns is false runList is a whitelist, otherwise a blacklist
  if (!tree)
    return kFALSE;

  if (!h) {
    Printf("Error - ProcessV0Tree: No THnSparse!");
    return kFALSE;
  }

  //set up the PID response
  if (!fPIDResponse) {
    fPIDResponse=new AliPIDResponse;
    fPIDResponse->SetOADBPath("$ALICE_PHYSICS/OADB");
    delete fESD;
    fESD = new AliESDEvent;
    fESD->CreateStdContent();
  }

  TClonesArray *arrTracks=(TClonesArray*)fESD->GetList()->FindObject("Tracks");

  // Tree name is important to assign PID
  TString treeName = tree->GetName();
  Bool_t isTreeK0 = treeName.CompareTo("treeK0", TString::kIgnoreCase) == 0;
  Bool_t isTreeGamma = treeName.CompareTo("treeGamma", TString::kIgnoreCase) == 0;
  Bool_t isTreeALambda = treeName.CompareTo("treeALambda", TString::kIgnoreCase) == 0;
  Bool_t isTreeLambda = treeName.CompareTo("treeLambda", TString::kIgnoreCase) == 0;

  if (!isTreeK0 && !isTreeGamma && !isTreeALambda && !isTreeLambda) {
    Printf("Error - ProcessV0Tree: Unknown tree type \"%s\"!", treeName.Data());
    return kFALSE;
  }

  const Int_t numCases = 5;
  Double_t tpcQA[numCases];
  Double_t tofQA[numCases];
  for (Int_t i = 0; i < numCases; i++) {
    tpcQA[i] = 0.;
    tofQA[i] = 0.;
  }

  const Int_t nTotESDTracks = -1;

  Int_t particleID = -1;
  AliPID::EParticleType alicePID=AliPID::kUnknown;

  Long64_t nTreeEntries = tree->GetEntriesFast();
  Double_t tpcsignal = -1.;
  Double_t tanTheta = -999.;
  Double_t precin = -1.;
  AliESDtrack* trk = 0x0;
  AliESDtrack* trk0 = 0x0;
  AliESDtrack* trk1 = 0x0;
  Int_t runNumber = -1;
  Int_t runNumberFirst = -1;
  Int_t ntracks   = -1;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("runNumber", 1);
  tree->SetBranchStatus("ntracks", 1);
  tree->SetBranchStatus("track0.*", 1);
  tree->SetBranchStatus("track1.*", 1);
  tree->SetBranchAddress("runNumber",&runNumber);
  tree->SetBranchAddress("ntracks"  ,&ntracks  );
  tree->SetBranchAddress("track0."  ,&trk0     );
  tree->SetBranchAddress("track1."  ,&trk1     );

  // Should be faster, but does not work for unknown reason
  //tree->SetBranchAddress(Form("track%d.fTPCsignal", track), &tpcsignal);
  //tree->SetBranchAddress(Form("track%d.fIp.fP[3]", track), &tanTheta);
  //tree->SetBranchAddress(Form("track%d.fIp.P()", track), &precin);

  // ===| set up run exclusion map |============================================
  TExMap runMap;
  if (!runList.IsNull()) {
    TObjArray *arr = runList.Tokenize(",; ");
    for (Int_t irun=0; irun<arr->GetEntriesFast(); ++irun) {
      const TString &runStr=((TObjString*)arr->At(irun))->String();
      if (!runStr.IsDigit()) continue;
      runMap.Add(runStr.Atoi(), 1);
      cout << "run" << runStr << " added" << endl;
    }
    delete arr;
  }

  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);
    // skip runs
    if (runMap.GetSize() && !(runMap.GetValue(runNumber)^excludeRuns)) continue;
    // set dummy esd track multiplicity
    arrTracks->Clear();
    arrTracks->ConstructedAt(ntracks);

    // initialise PID response
    // assume that we can load from the OADB only for the first run number
    // which means that the trees are merged per run or per period

    Bool_t init=kFALSE;
    if (runNumberFirst<0) {
      runNumberFirst=runNumber;
      init=kTRUE;
    }
    fPIDResponse->InitialiseEvent(fESD, recoPass, runNumberFirst);
    //if (init) {
      //fPIDResponse->GetTPCResponse().GetTransferFunctionParameters().Print();
    //}

    for (Int_t track = 0; track < 2; track++) {
      trk=trk0;
      if (track==1) {
        trk=trk1;
      }

      if (!trk)
        continue;
      //tpcsignal = trk->fTPCsignal;
      //tanTheta = trk->fIp->fP[3];
      //precin = trk->fIp->P();
      tpcsignal = trk->GetTPCsignal();
      tanTheta = trk->GetInnerParam()->GetParameter()[3];
      precin = trk->GetInnerParam()->P();

      //if (TMath::Abs(tanTheta-.7)>0.2) continue;
      // track == 0 is the positive particle, track == 1 the negative.
      // Combine this with the tree name to determine the particle ID.

      if (isTreeK0) {
        particleID = 2; //pion
        alicePID=AliPID::kPion;
      } else if (isTreeGamma) {
        particleID = 1; // electron
        alicePID=AliPID::kElectron;
      } else if (isTreeLambda) {
        if (track == 0) {
          particleID = 3; // proton
          alicePID=AliPID::kProton;
        } else {
          particleID = 2; // pion
          alicePID=AliPID::kPion;
        }
      }
      else if (isTreeALambda) {
        if (track == 0) {
          particleID = 2; // pion
          alicePID=AliPID::kPion;
        } else {
          particleID = 3; // proton
          alicePID=AliPID::kProton;
        }
      }

      // make signal correction for transfer function and BB shape
      // TODO: Perhaps this should be steerable
      //const Double_t corrTransferBB = fPIDResponse->GetTPCResponse().GetCombinedTransferFunctionBBCorrection(trk, alicePID);
      //const Double_t corrTransferBB = fPIDResponse->GetTPCResponse().GetTransferFunctionCorrection(trk);
      //printf("signal: %.2f/%.2f = %.2f\n", tpcsignal, corrTransferBB, tpcsignal/corrTransferBB);
      //tpcsignal/=corrTransferBB;

      Double_t processedTPCsignal[5] = { tpcsignal, tpcsignal, tpcsignal, tpcsignal, tpcsignal };
      /*
      if (fCorrectdEdxEtaDependence && fNumEtaCorrReqErrorsIssued < 23 && !etaCorrAvail) {
        AliError("TPC eta correction requested, but was not activated in PID response (most likely not available)!");
        fNumEtaCorrReqErrorsIssued++;
        if (fNumEtaCorrReqErrorsIssued == 23)
          AliError("Ignoring this error from now on!");
      }

      if (fCorrectdEdxMultiplicityDependence && fNumMultCorrReqErrorsIssued < 23 && !multCorrAvail) {
        AliError("TPC multiplicity correction requested, but was not activated in PID response (most likely not available)!");
        fNumMultCorrReqErrorsIssued++;
        if (fNumMultCorrReqErrorsIssued == 23)
          AliError("Ignoring this error from now on!");
      }

      if (fCorrectdEdxEtaDependence && etaCorrAvail && fCorrectdEdxMultiplicityDependence && multCorrAvail) {
        processedTPCsignal[0] = fPIDResponse->GetTPCResponse().GetEtaAndMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kElectron,
                                                                                                      AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[1] = fPIDResponse->GetTPCResponse().GetEtaAndMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kPion,
                                                                                                      AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[2] = fPIDResponse->GetTPCResponse().GetEtaAndMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kKaon,
                                                                                                      AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[3] = fPIDResponse->GetTPCResponse().GetEtaAndMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kProton,
                                                                                                      AliTPCPIDResponse::kdEdxDefault);
      }
      else if (fCorrectdEdxEtaDependence && etaCorrAvail) {
        processedTPCsignal[0] = fPIDResponse->GetTPCResponse().GetEtaCorrectedTrackdEdx(trackESD, AliPID::kElectron,
                                                                                        AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[1] = fPIDResponse->GetTPCResponse().GetEtaCorrectedTrackdEdx(trackESD, AliPID::kPion,
                                                                                        AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[2] = fPIDResponse->GetTPCResponse().GetEtaCorrectedTrackdEdx(trackESD, AliPID::kKaon,
                                                                                        AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[3] = fPIDResponse->GetTPCResponse().GetEtaCorrectedTrackdEdx(trackESD, AliPID::kProton,
                                                                                        AliTPCPIDResponse::kdEdxDefault);
      }
      else if (fCorrectdEdxMultiplicityDependence && multCorrAvail) {
        processedTPCsignal[0] = fPIDResponse->GetTPCResponse().GetMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kElectron,
                                                                                                AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[1] = fPIDResponse->GetTPCResponse().GetMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kPion,
                                                                                                AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[2] = fPIDResponse->GetTPCResponse().GetMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kKaon,
                                                                                                AliTPCPIDResponse::kdEdxDefault);
        processedTPCsignal[3] = fPIDResponse->GetTPCResponse().GetMultiplicityCorrectedTrackdEdx(trackESD, AliPID::kProton,
                                                                                                AliTPCPIDResponse::kdEdxDefault);
      }
      */

      // id 5 is just again Kaons in restricted eta range
      processedTPCsignal[4] = processedTPCsignal[2];

      for(Int_t iPart = 0; iPart < numCases; iPart++) {
        // Only accept "Kaons" within |eta| < 0.2 for index 4 (eta ~ tanTheta in this eta range)
        if (iPart == 4 && abs(tanTheta) > 0.2)
          continue;

        Double_t vecHistQA[7] = {precin, processedTPCsignal[iPart], (Double_t)particleID, (Double_t)iPart, tpcQA[iPart], tofQA[iPart],
                                (Double_t)nTotESDTracks};
        //cout << vecHistQA[0] << "\t" <<  vecHistQA[1] << "\t" <<  vecHistQA[2] << "\t" <<  vecHistQA[3] << "\t" <<  vecHistQA[4] << "\t" <<  vecHistQA[5] <<  endl;                      
        h->Fill(vecHistQA);
      }
    }
  }

  return kTRUE;
}


//________________________________________________________________________
THnSparseF* AliTPCcalibResidualPID::ProcessV0TreeFile(TString filePathName, const Int_t recoPass/*=4*/,
                                                      const TString runList/*=""*/, const Bool_t excludeRuns/*=kFALSE*/)
{
  // Open file and fill tracks from trees into a new THnSparse. Finally, returns it.

  TFile* f = TFile::Open(filePathName.Data(), "READ");
  if (!f) {
    Printf("Error - ProcessV0TreeFile: Cannot open file \"%s\"!", filePathName.Data());
    return 0x0;
  }

  THnSparseF* h = InitialisePIDQAHist("fHistPidQA","PID QA");

  TTree* tree = 0x0;
  const Int_t numTrees = 4;
  const TString treeNames[numTrees] = { "treeK0", "treeGamma", "treeLambda", "treeALambda" };
  for (Int_t i = 0; i < numTrees; i++) {
    Printf("Processing tree \"%s\"....", treeNames[i].Data());
    tree = (TTree*)f->Get(treeNames[i].Data());
    if (tree)
      tree->SetName(treeNames[i].Data());
    ProcessV0Tree(tree, h, recoPass, runList, excludeRuns);
  }

  f->Close();
  delete f;

  return h;
}


//________________________________________________________________________
Int_t  AliTPCcalibResidualPID::CompareFloat(Float_t f1, Float_t f2) const
{
    //compares if the Float_t f1 is equal with f2 and returns 1 if true and 0 if false
    Float_t precision = 0.00001;
    if (((f1 - precision) < f2) &&
        ((f1 + precision) > f2))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


//________________________________________________________________________
void AliTPCcalibResidualPID::Terminate(const Option_t *)
{


}



//________________________________________________________________________
void AliTPCcalibResidualPID::BinLogAxis(THnSparseF *h, Int_t axisNumber)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}

void AliTPCcalibResidualPID::CreatePlotWithOwnParameters(THnSparseF * histPidQA, const Bool_t useV0s, const Char_t * type, const Char_t * system, const Double_t* parameters, AliTPCcalibResidualPID::FitType fitType, Float_t from, Float_t to) {

  TObjArray * arr = 0x0;
  
  Bool_t isMC = kFALSE;
  if (strcmp(type, "MC") == 0) {
    arr = GetResidualGraphsMC(histPidQA, system);
    isMC = kTRUE;
  }
  else if (strcmp(type, "DATA") == 0) {
    arr = GetResidualGraphs(histPidQA, system, useV0s);
  }
  else {
    Printf("ERROR - ExtractResidualPID: Unknown type \"%s\" - must be \"MC\" or \"DATA\"!", type);
    
    return;
  }
  
  TF1* func = SetUpFitFunction(parameters, fitType, from, to, kFALSE, isMC, 0x0);
  
  TGraphErrors * graphAll = (TGraphErrors *) arr->FindObject("beamDataPoints");
  
  TCanvas* canvDelta_2 = CreateResidualCanvas(graphAll, func); 
  TCanvas* canvDelta_1 = CreateBBCanvas(arr, isMC, func); 
  return;
}

//________________________________________________________________________
Double_t* AliTPCcalibResidualPID::ExtractResidualPID(THnSparseF * histPidQA, const Bool_t useV0s, const Char_t * outFile,
                                                     const Char_t * type, const Char_t * period, const Char_t * pass,
                                                     const Char_t * system, const Double_t * initialParameters,
                                                     const Char_t *dedxtype,
                                                     AliTPCcalibResidualPID::FitType fitType) {
  //
  // (1.) obtain residual graphs
  //
  TObjArray * arr = 0x0;
  
  Bool_t isMC = kFALSE;
  if (strcmp(type, "MC") == 0) {
    arr = GetResidualGraphsMC(histPidQA, system);
    isMC = kTRUE;
  }
  else if (strcmp(type, "DATA") == 0) {
    arr = GetResidualGraphs(histPidQA, system, useV0s);
  }
  else {
    Printf("ERROR - ExtractResidualPID: Unknown type \"%s\" - must be \"MC\" or \"DATA\"!", type);
    
    return 0x0;
  }
  
  Bool_t isPPb = kFALSE;
  if (strcmp(system, "PPB") == 0 || strcmp(system, "PBP") == 0) {
    Printf("p-Pb/Pb-p detected - Applying special handling!");
    isPPb = kTRUE;
  }
  //
  // (2.) get old-style Bethe-Bloch parameters
  //
  TF1* parametrisation = FitBB(arr, isMC, isPPb, useV0s, initialParameters, fitType);
  //
  // (3.) obtain response functions to OADB
  //
  TObjArray* arrResponse = GetResponseFunctions(parametrisation, arr, type, period, pass, system, dedxtype);
  //
  // (4.) write results to file and exit
  //
  
  if (arrResponse) {
    TFile outputFile(outFile,"RECREATE");
    arrResponse->Write();
    outputFile.Close();
  }
  
  return parametrisation->GetParameters();
}


//________________________________________________________________________
TObjArray * AliTPCcalibResidualPID::GetSeparation(THnSparseF * histPidQA, Int_t kParticle1, Int_t kParticle2) {

  //get THnSparse
  THnSparse * hist = histPidQA;


  Float_t nSigmaPlus1 = 0, nSigmaMinus1 = 0; //sigma of first particle in the TPC
  Float_t nSigmaPlus2 = 0, nSigmaMinus2 = 0; //sigma of second particle in the TPC
  
  if(kParticle1 == 0){ // electron
        nSigmaMinus1   =  -1.999;
    nSigmaPlus1   =  2.999;
  } else if(kParticle1 == 1){ //pion
        nSigmaMinus1   =  -1.999;
    nSigmaPlus1   =  2.999;
  } else if(kParticle1 == 2){ //kaons
        nSigmaMinus1   =  -2.999;
    nSigmaPlus1   =    2.999;
  } else if(kParticle1 == 3){ //protons
        nSigmaMinus1   =  -2.999;
    nSigmaPlus1   =    2.999;
  }
  
  //
  // 1. select and fit first particle
  //
  hist->GetAxis(3)->SetRangeUser(kParticle1,kParticle1);
  hist->GetAxis(4)->SetRangeUser(nSigmaMinus1,nSigmaPlus1); 
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  TH2D * histPart1 = (TH2D*) hist->Projection(1,0)->Clone("histPart1");
  histPart1->SetTitle(";p /(GeV/c) ; TPCSignal /(a.u.)");
  histPart1->RebinX(4);

  TObjArray arr;
  histPart1->FitSlicesY(0,0,-1,10,"QNR",&arr);
  TH1D * PartMean1 = (TH1D *) arr.At(1)->Clone("PartMean1");
  TH1D * PartSigma1 = (TH1D *) arr.At(2)->Clone("PartSigma1");
  PartMean1->SetTitle(";p /(GeV/c) ; Mean (Gauss)");
  PartSigma1->SetTitle(";p /(GeV/c) ; Sigma (Gauss)");

  histPart1->SetMarkerStyle(22);
  PartMean1->SetMarkerStyle(21);
  hist->GetAxis(4)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES
  
  if(kParticle2==0){ // electron
        nSigmaMinus2   =  -1.999;
    nSigmaPlus2   =    2.999;
  } else if(kParticle2 == 1){ //pion
        nSigmaMinus2   =  -1.999;
    nSigmaPlus2   =    2.999;
  } else if(kParticle2 == 2){ //kaons
        nSigmaMinus2   =  -2.999;
    nSigmaPlus2   =    2.999;
  } else if(kParticle2 == 3){ //protons
        nSigmaMinus2   =  -2.999;
    nSigmaPlus2   =    2.999;
  }

  //
  // 2. select and fit second particle
  //
  hist->GetAxis(3)->SetRangeUser(kParticle2,kParticle2);
  hist->GetAxis(4)->SetRangeUser(nSigmaMinus2,nSigmaPlus2); 
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  TH2D * histPart2 = (TH2D*)hist->Projection(1,0)->Clone("histPart2");
  histPart2->RebinX(4);
  histPart2->FitSlicesY(0,0,-1,10,"QNR",&arr);
  TH1D * PartMean2 =  (TH1D *) arr.At(1)->Clone("PartMean2");
  TH1D * PartSigma2 =  (TH1D *) arr.At(2)->Clone("PartSigma2");
  PartMean2->SetTitle(";p /(GeV/c) ; Mean (Gauss)");
  PartSigma2->SetTitle(";p /(GeV/c) ; Sigma (Gauss)");
  PartMean2->SetMarkerStyle(20);
  hist->GetAxis(4)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES

  //
  // 3. separation
  //

        TH1F *fHistSeparation=(TH1F*)PartMean1->Clone("fHistSeparation"); //to get same binning
        fHistSeparation->SetMarkerStyle(22);
        const Int_t Nbins = PartMean1->GetNbinsX();
        fHistSeparation->SetTitle(";p /(GeV/c) ; Separation");

        Float_t DeltaMean[Nbins] ;
        Float_t DeltaSigma[Nbins];

        for(Int_t i=0 ; i<Nbins; i++){
        
        DeltaMean[i] = TMath::Abs(PartMean1->GetBinContent(i) - PartMean2->GetBinContent(i));
        DeltaSigma[i] = TMath::Abs((PartSigma1->GetBinContent(i) + PartSigma2->GetBinContent(i)))/2.;

        if(!(TMath::Abs(DeltaSigma[i])<0.000001))fHistSeparation->SetBinContent(i,DeltaMean[i]/DeltaSigma[i]);
        }//for(Int_t i=0 ; i<Nbins ; i++)

        TObjArray *array = new TObjArray();
        array->Add(histPart1);
        array->Add(histPart2);
        array->Add(PartMean1);
        array->Add(PartSigma1);
        array->Add(PartMean2);
        array->Add(PartSigma2);
        array->Add(fHistSeparation);
        return array;

}


//________________________________________________________________________
TObjArray * AliTPCcalibResidualPID::GetResidualGraphs(THnSparseF * histPidQA, const Char_t * system, Bool_t useV0s) {
  //
  // Extracts the residual graphs from THnSparse created from data (NOT MC!)
  //
  
  const TString momTitle = "#it{p}_{TPC} (GeV/#it{c})";
  const TString dEdxTitle = "d#it{E}/d#it{x} (arb. unit)";
  
  Bool_t isPPb = kFALSE;
  if (strcmp(system, "PPB") == 0 || strcmp(system, "PBP") == 0) {
    Printf("p-Pb/Pb-p detected - Applying special handling!");
    isPPb = kTRUE;
  }
  
  // the following parameters have to be extracted from the data itself
  //
    
  Float_t nSigmaPionMinus   =  -3.999;
  Float_t nSigmaPionPlus    =   3.999;
  Float_t nSigmaKaonMinus   =  -4.199;
  Float_t nSigmaKaonPlus    =  4.7999;
  Float_t nSigmaElectronMinus = -1.999;
  Float_t nSigmaElectronPlus  =  4.999;
  
  Int_t cutForFitting = 10;
  Double_t heightFractionForFittingRange = 0.1;
  //
  THnSparse * hist = histPidQA;
  //
  TCanvas * canvasQAtpc = new TCanvas("canvasQAtpcResGraph","Control canvas for residual graphs (TPC)",100,10,1380,800);
  canvasQAtpc->Divide(2,2);
  
  TCanvas * canvasQAtof = new TCanvas("canvasQAtofResGraph","Control canvas for residual graphs (TOF)",100,10,1380,800);
  canvasQAtof->Divide(2,2);
  
  TCanvas * canvasQAv0 = new TCanvas("canvasQAv0ResGraph","Control canvas for residual graphs (V0)",100,10,1380,800);
  canvasQAv0->Divide(2,2);
  
  TCanvas * canvasQAv0plusTOF = new TCanvas("canvasQAv0plusTOFResGraph","Control canvas for residual graphs (V0+TOF)",100,10,1380,800);
  canvasQAv0plusTOF->Divide(2,2);
  
  
  TCanvas * canvasQAv0DeDxPurityEl = new TCanvas("canvasQAv0DeDxPurityEl","Control canvas for residual graphs (V0 el)",100,10,600,400);
  TCanvas * canvasQAv0DeDxPurityPi = new TCanvas("canvasQAv0DeDxPurityPi","Control canvas for residual graphs (V0 pi)",100,10,600,400);
  TCanvas * canvasQAv0DeDxPurityPr = new TCanvas("canvasQAv0DeDxPurityPr","Control canvas for residual graphs (V0 pr)",100,10,600,400);
  
  canvasQAv0DeDxPurityEl->SetLogx();
  canvasQAv0DeDxPurityPi->SetLogx();
  canvasQAv0DeDxPurityPr->SetLogx();
  
  canvasQAv0DeDxPurityEl->SetLogz();
  canvasQAv0DeDxPurityPi->SetLogz();
  canvasQAv0DeDxPurityPr->SetLogz();
  
  canvasQAv0DeDxPurityEl->SetTopMargin(0.03);
  canvasQAv0DeDxPurityPi->SetTopMargin(0.03);
  canvasQAv0DeDxPurityPr->SetTopMargin(0.03);
  
  for (Int_t i = 1; i <= 4; i++) {
    canvasQAtpc->GetPad(i)->SetGrid(1, 1);
    canvasQAtof->GetPad(i)->SetGrid(1, 1);
    canvasQAv0->GetPad(i)->SetGrid(1, 1);
    canvasQAv0plusTOF->GetPad(i)->SetGrid(1, 1);

    canvasQAtpc->GetPad(i)->SetLogz();
    canvasQAtof->GetPad(i)->SetLogz();
    canvasQAv0->GetPad(i)->SetLogz();
    canvasQAv0plusTOF->GetPad(i)->SetLogz();

    canvasQAtpc->GetPad(i)->SetLogx();
    canvasQAtof->GetPad(i)->SetLogx();
    canvasQAv0->GetPad(i)->SetLogx();
    canvasQAv0plusTOF->GetPad(i)->SetLogx();
  }
  //
  // 1. select and fit electrons
  //
  hist->GetAxis(2)->SetRangeUser(0,0); // Non-V0s
  hist->GetAxis(3)->SetRangeUser(0,0); // electrons
  hist->GetAxis(4)->SetRangeUser(nSigmaElectronMinus,nSigmaElectronPlus); // 3sigma in TPC 
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  TH2D * histElectron = hist->Projection(1,0);
  histElectron->SetName("histElectron");
  histElectron->RebinX(4);
  histElectron->GetXaxis()->SetRangeUser(0.2,3.0);
  TObjArray arr;
  FitSlicesY(histElectron, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * electronPoints = (TH1D *) arr.At(1)->Clone();
  //
  TGraphErrors * electronGraph = new TGraphErrors(electronPoints);
  electronGraph->SetName("electronGraph");
  for(Int_t ip = 0; ip < electronGraph->GetN(); ip ++) {
    Bool_t removePoint = electronGraph->GetY()[ip] < 10 || 
      electronGraph->GetEY()[ip]/electronGraph->GetY()[ip] > 0.05 || 
      electronGraph->GetX()[ip] > 2.0 || 
      electronGraph->GetX()[ip] < 0.5; 
    if (removePoint) {
      electronGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAtof->cd(1);
  histElectron->SetTitle("electrons");
  histElectron->GetYaxis()->SetRangeUser(50, 120);
  histElectron->GetXaxis()->SetTitle(momTitle.Data());
  histElectron->GetYaxis()->SetTitle(dEdxTitle.Data());
  histElectron->GetXaxis()->SetMoreLogLabels(kTRUE);
  histElectron->GetXaxis()->SetNoExponent(kTRUE);
  histElectron->Draw("colz");
  electronPoints->SetMarkerStyle(24);
  electronPoints->Draw("same");
  hist->GetAxis(4)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES
  //
  // 2. protons
  //
  hist->GetAxis(2)->SetRangeUser(0,0); // Non-V0s
  hist->GetAxis(3)->SetRangeUser(3,3); // protons
  //
  //hist->GetAxis(4)->SetRangeUser(nSigmaProtonMinus,nSigmaProtonPlus); // protons --> not reliable, use PATTERN RECOGNITION
  TH2D * histProtonTPC = hist->Projection(1,0);
  histProtonTPC->SetName("histProtonTPC");
  histProtonTPC->GetXaxis()->SetRangeUser(0.15,0.7);
  histProtonTPC->GetYaxis()->SetRangeUser(50, hist->GetAxis(1)->GetBinUpEdge(hist->GetAxis(1)->GetNbins()));
  
  // PATTERN RECOGNITION
  //OLD TF1 betaSq("betaSq","50./TMath::Power(x,1.3)",0.1,10); 
  // Add some additional Erf - cuts away kaons for pPb and PbPb and does not hurt in pp
  TF1 betaSq("betaSq","50./TMath::Power(x,1.3) + (1-TMath::Erf((x-0.2) / 0.075)) * 100",0.1,10);
  for(Int_t ix = 1; ix <= histProtonTPC->GetXaxis()->GetNbins(); ix++) {
     for(Int_t jy = 1; jy <= histProtonTPC->GetYaxis()->GetNbins(); jy++) {
       Float_t yPos = histProtonTPC->GetYaxis()->GetBinCenter(jy);
       Float_t xPos = histProtonTPC->GetXaxis()->GetBinCenter(ix);
       Int_t bin = histProtonTPC->GetBin(ix,jy);
       if (yPos < betaSq.Eval(xPos)) histProtonTPC->SetBinContent(bin,0);
     }
  }
  //
  //TODO Use likelihood option -> Inside fitting range ~no contamination, but little statistics at low momenta requires L
  FitSlicesY(histProtonTPC, heightFractionForFittingRange, cutForFitting, "QNRL", &arr);
  TH1D * protonPointsTPC = (TH1D *) arr.At(1)->Clone();
  //
  TGraphErrors * protonTpcGraph = new TGraphErrors(protonPointsTPC);
  protonTpcGraph->SetName("protonTpcGraph");
  for(Int_t ip = 0; ip < protonTpcGraph->GetN(); ip ++) {
    Bool_t removePoint = protonTpcGraph->GetY()[ip] < 10 || protonTpcGraph->GetEY()[ip]/protonTpcGraph->GetY()[ip] > 0.05 // TODO Larger tolerance, since this is only for the correction function, where 5% are still better than no data point
                        || protonTpcGraph->GetX()[ip] > (useV0s ? (isPPb ? 0.499 : 0.349) : 0.649); // If V0s are to be used, don't use TPC only protons to too high momenta; for pPb low statistics for the moment, therefore, go a bit further with TPC protons
                        //TODO NOW -> Changed range of TPC-signal -> Hopefully, sufficient statistics now for very low momenta!|| protonTpcGraph->GetX()[ip] < 0.19; // Added this cut because almost always bad statistics below 0.19
    if (removePoint) {
      protonTpcGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  TH2D * histProtonTOF = hist->Projection(1,0);
  histProtonTOF->SetName("histProtonTOF");
  histProtonTOF->GetXaxis()->SetRangeUser(0.6,2.5);
  
  // Same pattern recognition as before
  for(Int_t ix = 1; ix <= histProtonTOF->GetXaxis()->GetNbins(); ix++) {
    for(Int_t jy = 1; jy <= histProtonTOF->GetYaxis()->GetNbins(); jy++) {
      Float_t yPos = histProtonTOF->GetYaxis()->GetBinCenter(jy);
      Float_t xPos = histProtonTOF->GetXaxis()->GetBinCenter(ix);
      Int_t bin = histProtonTOF->GetBin(ix,jy);
      if (yPos < betaSq.Eval(xPos)) histProtonTOF->SetBinContent(bin,0);
    }
  }
  
  FitSlicesY(histProtonTOF, heightFractionForFittingRange, cutForFitting, "QNR", &arr);

  TH1D * protonPointsTOF =  (TH1D *)arr.At(1)->Clone();
  //
  TGraphErrors * protonTofGraph = new TGraphErrors(protonPointsTOF);
  protonTofGraph->SetName("protonTofGraph");
  for(Int_t ip = 0; ip < protonTofGraph->GetN(); ip ++) {
    // If V0s are to be used, TPC+TOF protons are only for reference and can be plotted to higher momenta
    Bool_t removePoint = protonTofGraph->GetY()[ip] < 10 || protonTofGraph->GetEY()[ip]/protonTofGraph->GetY()[ip] > 0.02 || protonTofGraph->GetX()[ip] > (useV0s ? 3 : 2) || protonTofGraph->GetX()[ip] < 0.65;
    if (removePoint) {
      protonTofGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAtpc->cd(2);
  histProtonTPC->SetTitle("protons");
  histProtonTPC->GetXaxis()->SetTitle(momTitle.Data());
  histProtonTPC->GetYaxis()->SetTitle(dEdxTitle.Data());
  histProtonTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  histProtonTPC->GetXaxis()->SetNoExponent(kTRUE);
  histProtonTPC->Draw("colz");
  betaSq.DrawCopy("same");
  protonPointsTPC->SetMarkerStyle(20);
  protonPointsTOF->SetMarkerStyle(24);
  protonPointsTOF->Draw("same");
  protonPointsTPC->Draw("same");
  //
  protonTpcGraph->SetMarkerStyle(26);
  protonTpcGraph->SetMarkerColor(kMagenta);
  protonTpcGraph->DrawClone("p");
  protonTofGraph->SetMarkerStyle(25);
  protonTofGraph->SetMarkerColor(kMagenta);
  protonTofGraph->DrawClone("p");
  
  canvasQAtof->cd(2);
  histProtonTOF->SetTitle("protons");
  histProtonTOF->GetYaxis()->SetRangeUser(30, 250);
  histProtonTOF->GetXaxis()->SetTitle(momTitle.Data());
  histProtonTOF->GetYaxis()->SetTitle(dEdxTitle.Data());
  histProtonTOF->GetXaxis()->SetMoreLogLabels(kTRUE);
  histProtonTOF->GetXaxis()->SetNoExponent(kTRUE);
  histProtonTOF->Draw("colz");
  betaSq.DrawCopy("same");
  protonPointsTOF->Draw("same");
  protonPointsTPC->Draw("same");
  //
  protonTpcGraph->DrawClone("p");
  protonTofGraph->DrawClone("p");
  
  hist->GetAxis(4)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES
  //
  // 3. pions
  //
  hist->GetAxis(2)->SetRangeUser(0,0); // Non-V0s
  hist->GetAxis(3)->SetRangeUser(1,1); // pions
  //
  Double_t lowerPionThreshold = 0.3;

  hist->GetAxis(4)->SetRangeUser(nSigmaPionMinus,nSigmaPionPlus); // pions
  TH2D * histPionTPC = hist->Projection(1,0);
  histPionTPC->SetName("histPionTPC");
  histPionTPC->GetXaxis()->SetRangeUser(0.15,0.7);
  FitSlicesY(histPionTPC, heightFractionForFittingRange, cutForFitting, "QNR", &arr);

  // In case of really bad splines comment last 5 lines and use the following 6 lines:
  //TH2D * histPionTPC = hist->Projection(1,0);
  //histPionTPC->SetName("histPionTPC");
  //histPionTPC->GetYaxis()->SetRangeUser(30,100);
  //histPionTPC->GetXaxis()->SetRangeUser(0.15,0.7);
  //FitSlicesY(histPionTPC, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  //lowerPionThreshold = 0.15;


  TH1D * pionPointsTPC =  (TH1D *) arr.At(1)->Clone();
  //
  TGraphErrors * pionTpcGraph = new TGraphErrors(pionPointsTPC);
  pionTpcGraph->SetName("pionTpcGraph");
  for(Int_t ip = 0; ip < pionTpcGraph->GetN(); ip ++) {
    Bool_t removePoint = pionTpcGraph->GetY()[ip] < 10 || pionTpcGraph->GetEY()[ip]/pionTpcGraph->GetY()[ip] > 0.02 || pionTpcGraph->GetX()[ip] > 0.5
                         || pionTpcGraph->GetX()[ip] < lowerPionThreshold; // Exclude contamination by electrons
    if (removePoint) {
      pionTpcGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  TH2D * histPionTOF = hist->Projection(1,0);
  histPionTOF->SetName("histPionTOF");
  histPionTOF->GetXaxis()->SetRangeUser(0.5,1.1);
  FitSlicesY(histPionTOF, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * pionPointsTOF = (TH1D *) arr.At(1)->Clone();
  TGraphErrors * pionTofGraph = new TGraphErrors(pionPointsTOF);
  pionTofGraph->SetName("pionTofGraph");
  for(Int_t ip = 0; ip < pionTofGraph->GetN(); ip ++) {
    Bool_t removePoint = pionTofGraph->GetY()[ip] < 10 || pionTofGraph->GetEY()[ip]/pionTofGraph->GetY()[ip] > 0.02 || pionTofGraph->GetX()[ip] > 1.1 || pionTofGraph->GetX()[ip] < 0.5; // Changed by Ben (was 0.35) to avoid problems with TOF efficiency/systematics
    if (removePoint) {
      pionTofGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  canvasQAtpc->cd(3);
  histPionTPC->GetYaxis()->SetRangeUser(30, 90);
  histPionTPC->SetTitle("pions");
  histPionTPC->GetXaxis()->SetTitle(momTitle.Data());
  histPionTPC->GetYaxis()->SetTitle(dEdxTitle.Data());
  histPionTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  histPionTPC->GetXaxis()->SetNoExponent(kTRUE);
  histPionTPC->Draw("colz");
  pionPointsTPC->SetMarkerStyle(20);
  pionPointsTOF->SetMarkerStyle(24);
  pionPointsTOF->Draw("same");
  pionPointsTPC->Draw("same");
  //
  pionTpcGraph->SetMarkerStyle(26);
  pionTpcGraph->SetMarkerColor(kMagenta);
  pionTpcGraph->DrawClone("p");
  pionTofGraph->SetMarkerStyle(25);
  pionTofGraph->SetMarkerColor(kMagenta);
  pionTofGraph->DrawClone("p");
  
  canvasQAtof->cd(3);
  histPionTOF->GetYaxis()->SetRangeUser(30, 80);
  histPionTOF->GetXaxis()->SetTitle(momTitle.Data());
  histPionTOF->GetYaxis()->SetTitle(dEdxTitle.Data());
  histPionTOF->GetXaxis()->SetMoreLogLabels(kTRUE);
  histPionTOF->GetXaxis()->SetNoExponent(kTRUE);
  histPionTOF->Draw("colz");
  histPionTOF->SetTitle("pions");
  pionPointsTOF->Draw("same");
  pionPointsTPC->Draw("same");
  //
  pionTpcGraph->DrawClone("p");
  pionTofGraph->DrawClone("p");

  
  hist->GetAxis(4)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES
  //
  // 4. kaons
  //
  hist->GetAxis(2)->SetRangeUser(0,0); // Non-V0s
  hist->GetAxis(3)->SetRangeUser(2,2); // kaons
  //
  hist->GetAxis(4)->SetRangeUser(nSigmaKaonMinus,nSigmaKaonPlus); // kaons
  TH2D * histKaonTPC = hist->Projection(1,0);
  histKaonTPC->SetName("histKaonTPC");
  histKaonTPC->GetXaxis()->SetRangeUser(0.12,0.6999);
  FitSlicesY(histKaonTPC, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * kaonPointsTPC = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * kaonTpcGraph = new TGraphErrors(kaonPointsTPC);
  kaonTpcGraph->SetName("kaonTpcGraph");
  for(Int_t ip = 0; ip < kaonTpcGraph->GetN(); ip ++) {
    Bool_t removePoint = kaonTpcGraph->GetY()[ip] < 10 || kaonTpcGraph->GetEY()[ip]/kaonTpcGraph->GetY()[ip] > 0.02
                         || kaonTpcGraph->GetX()[ip] < 0.12 || kaonTpcGraph->GetX()[ip] > 0.319; //TODO BEN was > 0.399
    if (removePoint) {
      kaonTpcGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // sigma in TOF
  
  /*
  hist->GetAxis(5)->SetRangeUser(-1.999,1.999); // sigma in TOF
  
  if (hist->GetAxis(3)->GetNbins() >= 5) {
    printf("Using restricted eta range for TPC+TOF Kaons!\n");
    hist->GetAxis(3)->SetRangeUser(4,4);
  }
  else {
    printf("WARNING: Data points for restricted eta range for TPC+TOF Kaons not available! Using full eta range (worse w.r.t. contamination)\n");
  }*/
  
  TH2D * histKaonTOF = hist->Projection(1,0);
  histKaonTOF->SetName("histKaonTOF");
  histKaonTOF->GetXaxis()->SetRangeUser(0.3,1.5);
  FitSlicesY(histKaonTOF, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * kaonPointsTOF = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * kaonTofGraph = new TGraphErrors(kaonPointsTOF);
  kaonTofGraph->SetName("kaonTofGraph");
  //
  for(Int_t ip = 0; ip < kaonTofGraph->GetN(); ip ++) {
    Bool_t removePoint = kaonTofGraph->GetY()[ip] < 10 || kaonTofGraph->GetEY()[ip]/kaonTofGraph->GetY()[ip] > 0.02 || kaonTofGraph->GetX()[ip] > 1.0
                         || kaonTofGraph->GetX()[ip] < 0.36; //TODO BEN NOW was 0.5
    if (removePoint) {
      kaonTofGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAtpc->cd(4);
  histKaonTPC->SetTitle("kaons");
  histKaonTPC->GetYaxis()->SetRangeUser(30, 500);
  histKaonTPC->GetXaxis()->SetTitle(momTitle.Data());
  histKaonTPC->GetYaxis()->SetTitle(dEdxTitle.Data());
  histKaonTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  histKaonTPC->GetXaxis()->SetNoExponent(kTRUE);
  histKaonTPC->Draw("colz");
  kaonPointsTPC->SetMarkerStyle(20);
  kaonPointsTOF->SetMarkerStyle(24);
  kaonPointsTOF->Draw("same");
  kaonPointsTPC->Draw("same");
  //
  kaonTpcGraph->SetMarkerStyle(26);
  kaonTpcGraph->SetMarkerColor(kMagenta);
  kaonTpcGraph->DrawClone("p");
  kaonTofGraph->SetMarkerStyle(25);
  kaonTofGraph->SetMarkerColor(kMagenta);
  kaonTofGraph->DrawClone("p");

  
  canvasQAtof->cd(4);
  histKaonTOF->GetYaxis()->SetRangeUser(30, 250);
  histKaonTOF->SetTitle("kaons");
  histKaonTOF->GetXaxis()->SetTitle(momTitle.Data());
  histKaonTOF->GetYaxis()->SetTitle(dEdxTitle.Data());
  histKaonTOF->GetXaxis()->SetMoreLogLabels(kTRUE);
  histKaonTOF->GetXaxis()->SetNoExponent(kTRUE);
  histKaonTOF->Draw("colz");
  kaonPointsTOF->Draw("same");
  kaonPointsTPC->Draw("same");
  kaonTpcGraph->DrawClone("p");
  kaonTofGraph->DrawClone("p");
  
  hist->GetAxis(4)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES
  
  
  //
  // 5. V0 electrons
  //
  hist->GetAxis(2)->SetRangeUser(1,1); // V0 electrons
  hist->GetAxis(3)->SetRangeUser(0,0); // electrons
  
  //
  TH2D * histElectronV0 = hist->Projection(1,0);
  histElectronV0->SetName("histElectronV0");
  //histElectronV0->RebinX(2);
  histElectronV0->GetXaxis()->SetRangeUser(0.15,5.0);
  FitSlicesY(histElectronV0, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * electronPointsV0 = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * electronV0Graph = new TGraphErrors(electronPointsV0);
  electronV0Graph->SetName("electronV0Graph");
  for(Int_t ip = 0; ip < electronV0Graph->GetN(); ip ++) {
    Bool_t removePoint = electronV0Graph->GetY()[ip] < 10 || electronV0Graph->GetEY()[ip]/electronV0Graph->GetY()[ip] > 0.02;

    if (removePoint) {
      electronV0Graph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAv0->cd(1);
  histElectronV0->SetTitle("V0 electrons");
  histElectronV0->GetYaxis()->SetRangeUser(50, 120);
  histElectronV0->GetXaxis()->SetTitle(momTitle.Data());
  histElectronV0->GetYaxis()->SetTitle(dEdxTitle.Data());
  histElectronV0->GetXaxis()->SetMoreLogLabels(kTRUE);
  histElectronV0->GetXaxis()->SetNoExponent(kTRUE);
  histElectronV0->SetStats(0);
  histElectronV0->Draw("colz");
  electronPointsV0->SetMarkerStyle(34);
  electronPointsV0->Draw("same");
  electronGraph->SetMarkerStyle(25);
  electronGraph->SetMarkerColor(kMagenta);
  electronGraph->DrawClone("p");
  
  canvasQAv0DeDxPurityEl->cd();
  gStyle->SetOptTitle(0);
  TH1D* hElNew = (TH1D*)histElectronV0->DrawClone("colz");
  hElNew->GetXaxis()->SetTitleOffset(1.0);
  hElNew->SetTitle("");
  electronPointsV0->DrawClone("same");
  gStyle->SetOptTitle(1);



  canvasQAtof->cd(1);
  electronV0Graph->SetMarkerStyle(24);
  electronV0Graph->SetMarkerColor(kMagenta);
  electronV0Graph->DrawClone("p");
  
  //
  // 6. V0 pions
  //
  hist->GetAxis(2)->SetRangeUser(2,2); // V0 pions
  hist->GetAxis(3)->SetRangeUser(1,1); // pions
  //
  TH2D * histPionV0 = hist->Projection(1,0);
  histPionV0->SetName("histPionV0");
  histPionV0->GetXaxis()->SetRangeUser(0.15,10.);
  FitSlicesY(histPionV0, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * pionPointsV0 = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * pionV0Graph = new TGraphErrors(pionPointsV0);
  pionV0Graph->SetName("pionV0Graph");
  for(Int_t ip = 0; ip < pionV0Graph->GetN(); ip ++) {
    Bool_t removePoint = pionV0Graph->GetY()[ip] < 10 || pionV0Graph->GetEY()[ip]/pionV0Graph->GetY()[ip] > 0.02;

    if (removePoint) {
      pionV0Graph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAv0->cd(3);
  histPionV0->SetTitle("V0 pions");
  histPionV0->GetYaxis()->SetRangeUser(30, 100);
  histPionV0->GetXaxis()->SetTitle(momTitle.Data());
  histPionV0->GetYaxis()->SetTitle(dEdxTitle.Data());
  histPionV0->GetXaxis()->SetMoreLogLabels(kTRUE);
  histPionV0->GetXaxis()->SetNoExponent(kTRUE);
  histPionV0->Draw("colz");
  pionPointsV0->SetMarkerStyle(34);
  pionPointsV0->Draw("same");
  pionTofGraph->DrawClone("p");
  pionTpcGraph->DrawClone("p");
  
  canvasQAv0DeDxPurityPi->cd();
  gStyle->SetOptTitle(0);
  TH1D* hPiNew = (TH1D*)histPionV0->DrawClone("colz");
  hPiNew->GetXaxis()->SetTitleOffset(1.0);
  hPiNew->SetTitle("");
  pionPointsV0->DrawClone("same");
  gStyle->SetOptTitle(1);
  
  canvasQAtof->cd(3);
  pionV0Graph->SetMarkerStyle(24);
  pionV0Graph->SetMarkerColor(kMagenta);
  pionV0Graph->DrawClone("p");

  canvasQAtpc->cd(3);
  pionV0Graph->DrawClone("p");
  
  //
  // 6. V0 protons
  //
  hist->GetAxis(2)->SetRangeUser(3,3); // V0 protons
  hist->GetAxis(3)->SetRangeUser(3,3); // protons
  //
  TH2D * histProtonV0 = hist->Projection(1,0);
  histProtonV0->SetName("histProtonV0");
  histProtonV0->GetXaxis()->SetRangeUser(0.15,10.);
  
  // Remove contamination due to el and pions and low momenta because there the statistics
  // for the V0 protons goes down and becomes comparable with the contamination
  // -> This spoils the fit, if the contamination is not removed
  
  // PATTERN RECOGNITION
  Int_t correctionUpToBin = histProtonV0->GetXaxis()->FindBin(0.6);
  
  for(Int_t ix = 1; ix <= correctionUpToBin; ix++) {
     for(Int_t jy = 1; jy <= histProtonV0->GetYaxis()->GetNbins(); jy++) {
       Float_t yPos = histProtonV0->GetYaxis()->GetBinCenter(jy);
       Float_t xPos = histProtonV0->GetXaxis()->GetBinCenter(ix);
       Int_t bin = histProtonV0->GetBin(ix,jy);
       if (yPos < betaSq.Eval(xPos)) histProtonV0->SetBinContent(bin,0);
     }
  }
  
  /*
  Int_t correctionUpToBin = histProtonV0->GetXaxis()->FindBin(0.4);
  Double_t threshold = 120.;
  
  for(Int_t ix = 1; ix <= correctionUpToBin; ix++) {
     for(Int_t jy = 1; jy <= histProtonV0->GetYaxis()->GetNbins(); jy++) {
       Float_t yPos = histProtonV0->GetYaxis()->GetBinCenter(jy);
       Int_t bin = histProtonV0->GetBin(ix,jy);
       if (yPos < threshold) histProtonV0->SetBinContent(bin,0);
       else break;
     }
  }
  */
  
  FitSlicesY(histProtonV0, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * protonPointsV0 = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * protonV0Graph = new TGraphErrors(protonPointsV0);
  protonV0Graph->SetName("protonV0Graph");
  for(Int_t ip = 0; ip < protonV0Graph->GetN(); ip ++) {
    Bool_t removePoint = protonV0Graph->GetY()[ip] < 10 || protonV0Graph->GetEY()[ip]/protonV0Graph->GetY()[ip] > 0.02;
                         
    if (removePoint) {
      protonV0Graph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAv0->cd(2);
  histProtonV0->SetTitle("V0 protons");
  histProtonV0->GetXaxis()->SetTitle(momTitle.Data());
  histProtonV0->GetYaxis()->SetTitle(dEdxTitle.Data());
  histProtonV0->GetXaxis()->SetMoreLogLabels(kTRUE);
  histProtonV0->GetXaxis()->SetNoExponent(kTRUE);
  histProtonV0->Draw("colz");
  protonPointsV0->SetMarkerStyle(34);
  protonPointsV0->Draw("same");
  protonTofGraph->DrawClone("p");
  protonTpcGraph->DrawClone("p");
  
  canvasQAv0DeDxPurityPr->cd();
  gStyle->SetOptTitle(0);
  TH1D* hPrNew = (TH1D*)histProtonV0->DrawClone("colz");
  hPrNew->GetXaxis()->SetTitleOffset(1.0);
  hPrNew->SetTitle("");
  protonPointsV0->DrawClone("same");
  gStyle->SetOptTitle(1);

  canvasQAtof->cd(2);
  protonV0Graph->SetMarkerStyle(24);
  protonV0Graph->SetMarkerColor(kMagenta);
  protonV0Graph->DrawClone("p");

  canvasQAtpc->cd(2);
  protonV0Graph->DrawClone("p");
  
  
  hist->GetAxis(2)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(3)->SetRange(0,-1); // RESET RANGES
  
  
  //
  // 5. V0 electrons + TOF
  //
  hist->GetAxis(2)->SetRangeUser(1,1); // V0 electrons
  hist->GetAxis(3)->SetRangeUser(0,0); // electrons
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  //
  TH2D * histElectronV0plusTOF = hist->Projection(1,0);
  histElectronV0plusTOF->SetName("histElectronV0plusTOF");
  histElectronV0plusTOF->RebinX(2);
  histElectronV0plusTOF->GetXaxis()->SetRangeUser(0.2,5.0);
  FitSlicesY(histElectronV0plusTOF, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * electronPointsV0plusTOF = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * electronV0plusTOFGraph = new TGraphErrors(electronPointsV0plusTOF);
  electronV0plusTOFGraph->SetName("electronV0plusTOFGraph");
  for(Int_t ip = 0; ip < electronV0plusTOFGraph->GetN(); ip ++) {
    Bool_t removePoint = electronV0plusTOFGraph->GetY()[ip] < 10 || electronV0plusTOFGraph->GetEY()[ip]/electronV0plusTOFGraph->GetY()[ip] > 0.02;

    if (removePoint) {
      electronV0plusTOFGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAv0plusTOF->cd(1);
  histElectronV0plusTOF->SetTitle("V0+TOF electrons");
  histElectronV0plusTOF->GetYaxis()->SetRangeUser(50, 120);
  histElectronV0plusTOF->GetXaxis()->SetTitle(momTitle.Data());
  histElectronV0plusTOF->GetYaxis()->SetTitle(dEdxTitle.Data());
  histElectronV0plusTOF->GetXaxis()->SetMoreLogLabels(kTRUE);
  histElectronV0plusTOF->GetXaxis()->SetNoExponent(kTRUE);
  histElectronV0plusTOF->Draw("colz");
  electronPointsV0plusTOF->SetMarkerStyle(29);
  electronPointsV0plusTOF->Draw("same");
  electronGraph->DrawClone("p");
  electronV0Graph->DrawClone("p");

  canvasQAv0->cd(1);
  electronV0plusTOFGraph->SetMarkerStyle(30);
  electronV0plusTOFGraph->SetMarkerColor(kMagenta);
  electronV0plusTOFGraph->DrawClone("p");

  canvasQAtof->cd(1);
  electronV0plusTOFGraph->DrawClone("p");
  
  //
  // 6. V0 pions + TOF
  //
  hist->GetAxis(2)->SetRangeUser(2,2); // V0 pions
  hist->GetAxis(3)->SetRangeUser(1,1); // pions
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  //
  TH2D * histPionV0plusTOF = hist->Projection(1,0);
  histPionV0plusTOF->SetName("histPionV0plusTOF");
  histPionV0plusTOF->GetXaxis()->SetRangeUser(0.15,10.);
  FitSlicesY(histPionV0plusTOF, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * pionPointsV0plusTOF = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * pionV0plusTOFGraph = new TGraphErrors(pionPointsV0plusTOF);
  pionV0plusTOFGraph->SetName("pionV0plusTOFGraph");
  for(Int_t ip = 0; ip < pionV0plusTOFGraph->GetN(); ip ++) {
    Bool_t removePoint = pionV0plusTOFGraph->GetY()[ip] < 10 || pionV0plusTOFGraph->GetEY()[ip]/pionV0plusTOFGraph->GetY()[ip] > 0.02;

    if (removePoint) {
      pionV0plusTOFGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAv0plusTOF->cd(3);
  histPionV0plusTOF->SetTitle("V0+TOF pions");
  histPionV0plusTOF->GetYaxis()->SetRangeUser(30, 100);
  histPionV0plusTOF->GetXaxis()->SetTitle(momTitle.Data());
  histPionV0plusTOF->GetYaxis()->SetTitle(dEdxTitle.Data());
  histPionV0plusTOF->GetXaxis()->SetMoreLogLabels(kTRUE);
  histPionV0plusTOF->GetXaxis()->SetNoExponent(kTRUE);
  histPionV0plusTOF->Draw("colz");
  pionPointsV0plusTOF->SetMarkerStyle(29);
  pionPointsV0plusTOF->Draw("same");
  pionV0Graph->DrawClone("p");
  pionTofGraph->DrawClone("p");
  pionTpcGraph->DrawClone("p");
  
  canvasQAv0->cd(3);
  pionV0plusTOFGraph->SetMarkerStyle(30);
  pionV0plusTOFGraph->SetMarkerColor(kMagenta);
  pionV0plusTOFGraph->DrawClone("p");
  
  canvasQAtof->cd(3);
  pionV0plusTOFGraph->DrawClone("p");
  
  canvasQAtpc->cd(3);
  pionV0plusTOFGraph->DrawClone("p");
  
  //
  // 6. V0 protons + TOF
  //
  hist->GetAxis(2)->SetRangeUser(3,3); // V0 protons
  hist->GetAxis(3)->SetRangeUser(3,3); // protons
  hist->GetAxis(5)->SetRangeUser(-2.999,2.999); // 3sigma in TOF
  //
  TH2D * histProtonV0plusTOF = hist->Projection(1,0);
  histProtonV0plusTOF->SetName("histProtonV0plusTOF");
  histProtonV0plusTOF->GetXaxis()->SetRangeUser(0.15,10.);
  
  // Remove contamination due to el and pions and low momenta because there the statistics
  // for the V0 protons goes down and becomes comparable with the contamination
  // -> This spoils the fit, if the contamination is not removed
  
  // PATTERN RECOGNITION
  correctionUpToBin = histProtonV0plusTOF->GetXaxis()->FindBin(0.6);
  
  for(Int_t ix = 1; ix <= correctionUpToBin; ix++) {
     for(Int_t jy = 1; jy <= histProtonV0plusTOF->GetYaxis()->GetNbins(); jy++) {
       Float_t yPos = histProtonV0plusTOF->GetYaxis()->GetBinCenter(jy);
       Float_t xPos = histProtonV0plusTOF->GetXaxis()->GetBinCenter(ix);
       Int_t bin = histProtonV0plusTOF->GetBin(ix,jy);
       if (yPos < betaSq.Eval(xPos)) histProtonV0plusTOF->SetBinContent(bin,0);
     }
  }
  
  FitSlicesY(histProtonV0plusTOF, heightFractionForFittingRange, cutForFitting, "QNR", &arr);
  TH1D * protonPointsV0plusTOF = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * protonV0plusTOFGraph = new TGraphErrors(protonPointsV0plusTOF);
  protonV0plusTOFGraph->SetName("protonV0plusTOFGraph");
  for(Int_t ip = 0; ip < protonV0plusTOFGraph->GetN(); ip ++) {
    Bool_t removePoint = protonV0plusTOFGraph->GetY()[ip] < 10 || protonV0plusTOFGraph->GetEY()[ip]/protonV0plusTOFGraph->GetY()[ip] > 0.02;
                         
    if (removePoint) {
      protonV0plusTOFGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAv0plusTOF->cd(2);
  histProtonV0plusTOF->SetTitle("V0+TOF protons");
  histProtonV0plusTOF->GetXaxis()->SetTitle(momTitle.Data());
  histProtonV0plusTOF->GetYaxis()->SetTitle(dEdxTitle.Data());
  histProtonV0plusTOF->GetXaxis()->SetMoreLogLabels(kTRUE);
  histProtonV0plusTOF->GetXaxis()->SetNoExponent(kTRUE);
  histProtonV0plusTOF->Draw("colz");
  protonPointsV0plusTOF->SetMarkerStyle(29);
  protonPointsV0plusTOF->Draw("same");
  protonV0Graph->DrawClone("p");
  protonTofGraph->DrawClone("p");
  protonTpcGraph->DrawClone("p");

  canvasQAv0->cd(2);
  protonV0plusTOFGraph->SetMarkerStyle(30);
  protonV0plusTOFGraph->SetMarkerColor(kMagenta);
  protonV0plusTOFGraph->DrawClone("p");
  
  canvasQAtof->cd(2);
  protonV0plusTOFGraph->DrawClone("p");

  canvasQAtpc->cd(2);
  protonV0plusTOFGraph->DrawClone("p");
  
  hist->GetAxis(2)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(3)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(5)->SetRange(0,-1); // RESET RANGES
  
  
  // Clone (V0, V0+TOF) graphs: Remove some points, which are expected to deviate from the Bethe-Bloch curve, from the original graphs, so that
  // these graphs can be used for the BB fit. The original graph will keep all point in order to determine the correction and response
  // functions
  TGraphErrors * electronV0GraphForBBfit = new TGraphErrors(*electronV0Graph);
  electronV0GraphForBBfit->SetName("electronV0GraphForBBfit");
  
  for(Int_t ip = 0; ip < electronV0GraphForBBfit->GetN(); ip ++) {
    Bool_t removePoint = electronV0GraphForBBfit->GetX()[ip] < (isPPb ? 1.2 : 0.4);// || electronV0GraphForBBfit->GetX()[ip] > 2.2;
    
    if (removePoint) {
      electronV0GraphForBBfit->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  TGraphErrors * pionV0GraphForBBfit = new TGraphErrors(*pionV0Graph);
  pionV0GraphForBBfit->SetName("pionV0GraphForBBfit");
  
  for(Int_t ip = 0; ip < pionV0GraphForBBfit->GetN(); ip ++) {
    Bool_t removePoint = pionV0GraphForBBfit->GetX()[ip] < (isPPb ? 1.0 : 0.76);// || pionV0GraphForBBfit->GetX()[ip] > 6.; //TODO NOW was < 0.4 before
    //Bool_t removePoint = pionV0GraphForBBfit->GetX()[ip] < 0.4 || pionV0GraphForBBfit->GetX()[ip] > 0.399;// In this case NOT used for BB fit
    if (removePoint) {
      pionV0GraphForBBfit->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  TGraphErrors * protonV0GraphForBBfit = new TGraphErrors(*protonV0Graph);
  protonV0GraphForBBfit->SetName("protonV0GraphForBBfit");
  
  for(Int_t ip = 0; ip < protonV0GraphForBBfit->GetN(); ip ++) {
    Bool_t removePoint = protonV0GraphForBBfit->GetX()[ip] < 0.9 || protonV0GraphForBBfit->GetX()[ip] > 5.0; //TODO NOW was 2.5 before
    //Bool_t removePoint = protonV0GraphForBBfit->GetX()[ip] < 0.9 || protonV0GraphForBBfit->GetX()[ip] > 0.599;// In this case NOT used for BB fit
    if (removePoint) {
      protonV0GraphForBBfit->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  
  // V0 plus TOF
  TGraphErrors * electronV0plusTOFGraphForBBfit = new TGraphErrors(*electronV0plusTOFGraph);
  electronV0plusTOFGraphForBBfit->SetName("electronV0plusTOFGraphForBBfit");
  
  for(Int_t ip = 0; ip < electronV0plusTOFGraphForBBfit->GetN(); ip ++) {
    Bool_t removePoint = electronV0plusTOFGraphForBBfit->GetX()[ip] < (isPPb ? 1.2 : 0.6);//0.4 || electronV0plusTOFGraphForBBfit->GetX()[ip] > 1.3799;
    if (removePoint) {
      electronV0plusTOFGraphForBBfit->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  TGraphErrors * pionV0plusTOFGraphForBBfit = new TGraphErrors(*pionV0plusTOFGraph);
  pionV0plusTOFGraphForBBfit->SetName("pionV0plusTOFGraphForBBfit");
  
  for(Int_t ip = 0; ip < pionV0plusTOFGraphForBBfit->GetN(); ip ++) {
    Bool_t removePoint = pionV0plusTOFGraphForBBfit->GetX()[ip] < 0.4;// || pionV0plusTOFGraphForBBfit->GetX()[ip] > 6.;
    if (removePoint) {
      pionV0plusTOFGraphForBBfit->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  TGraphErrors * protonV0plusTOFGraphForBBfit = new TGraphErrors(*protonV0plusTOFGraph);
  protonV0plusTOFGraphForBBfit->SetName("protonV0plusTOFGraphForBBfit");
  
  for(Int_t ip = 0; ip < protonV0plusTOFGraphForBBfit->GetN(); ip ++) {
    Bool_t removePoint = protonV0plusTOFGraphForBBfit->GetX()[ip] < 0.9 || protonV0plusTOFGraphForBBfit->GetX()[ip] > 2.5;
    if (removePoint) {
      protonV0plusTOFGraphForBBfit->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  
  //
  // rescale graphs to betaGamma
  //
  kaonTofGraph->SetMarkerStyle(21);
  kaonTpcGraph->SetMarkerStyle(22);
  for(Int_t ip = 0; ip < kaonTofGraph->GetN(); ip ++) {
    kaonTofGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
    kaonTofGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
  }
  for(Int_t ip = 0; ip < kaonTpcGraph->GetN(); ip ++) {
    kaonTpcGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
    kaonTpcGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
  }
  //
  electronGraph->SetMarkerStyle(22);
  electronV0Graph->SetMarkerStyle(34);
  electronV0GraphForBBfit->SetMarkerStyle(34);
  electronV0plusTOFGraph->SetMarkerStyle(30);
  electronV0plusTOFGraphForBBfit->SetMarkerStyle(30);
  electronGraph->SetMarkerColor(2);
  electronV0Graph->SetMarkerColor(2);
  electronV0GraphForBBfit->SetMarkerColor(2);
  electronV0plusTOFGraph->SetMarkerColor(2);
  electronV0plusTOFGraphForBBfit->SetMarkerColor(2);
  for(Int_t ip = 0; ip < electronGraph->GetN(); ip ++) {
    electronGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    electronGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  for(Int_t ip = 0; ip < electronV0Graph->GetN(); ip ++) {
    electronV0Graph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    electronV0Graph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  for(Int_t ip = 0; ip < electronV0GraphForBBfit->GetN(); ip ++) {
    electronV0GraphForBBfit->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    electronV0GraphForBBfit->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  for(Int_t ip = 0; ip < electronV0plusTOFGraph->GetN(); ip ++) {
    electronV0plusTOFGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    electronV0plusTOFGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  for(Int_t ip = 0; ip < electronV0plusTOFGraphForBBfit->GetN(); ip ++) {
    electronV0plusTOFGraphForBBfit->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    electronV0plusTOFGraphForBBfit->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  //
  pionTofGraph->SetMarkerStyle(21);
  pionTpcGraph->SetMarkerStyle(22);
  pionV0Graph->SetMarkerStyle(34);
  pionV0GraphForBBfit->SetMarkerStyle(34);
  pionV0plusTOFGraph->SetMarkerStyle(30);
  pionV0plusTOFGraphForBBfit->SetMarkerStyle(30);
  pionTofGraph->SetMarkerColor(3);
  pionTpcGraph->SetMarkerColor(3);
  pionV0Graph->SetMarkerColor(3);
  pionV0GraphForBBfit->SetMarkerColor(3);
  pionV0plusTOFGraph->SetMarkerColor(3);
  pionV0plusTOFGraphForBBfit->SetMarkerColor(3);
  for(Int_t ip = 0; ip < pionTofGraph->GetN(); ip ++) {
    pionTofGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionTofGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  for(Int_t ip = 0; ip < pionTpcGraph->GetN(); ip ++) {
    pionTpcGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionTpcGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  for(Int_t ip = 0; ip < pionV0Graph->GetN(); ip ++) {
    pionV0Graph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionV0Graph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  for(Int_t ip = 0; ip < pionV0GraphForBBfit->GetN(); ip ++) {
    pionV0GraphForBBfit->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionV0GraphForBBfit->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  for(Int_t ip = 0; ip < pionV0plusTOFGraph->GetN(); ip ++) {
    pionV0plusTOFGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionV0plusTOFGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  for(Int_t ip = 0; ip < pionV0plusTOFGraphForBBfit->GetN(); ip ++) {
    pionV0plusTOFGraphForBBfit->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionV0plusTOFGraphForBBfit->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  //
  protonTofGraph->SetMarkerStyle(21);
  protonTpcGraph->SetMarkerStyle(22);
  protonV0Graph->SetMarkerStyle(34);
  protonV0GraphForBBfit->SetMarkerStyle(34);
  protonV0plusTOFGraph->SetMarkerStyle(30);
  protonV0plusTOFGraphForBBfit->SetMarkerStyle(30);
  protonTofGraph->SetMarkerColor(4);
  protonTpcGraph->SetMarkerColor(4);
  protonV0Graph->SetMarkerColor(4);
  protonV0GraphForBBfit->SetMarkerColor(4);
  protonV0plusTOFGraph->SetMarkerColor(4);
  protonV0plusTOFGraphForBBfit->SetMarkerColor(4);
  for(Int_t ip = 0; ip < protonTofGraph->GetN(); ip ++) {
    protonTofGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonTofGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  for(Int_t ip = 0; ip < protonTpcGraph->GetN(); ip ++) {
    protonTpcGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonTpcGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  for(Int_t ip = 0; ip < protonV0Graph->GetN(); ip ++) {
    protonV0Graph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonV0Graph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  for(Int_t ip = 0; ip < protonV0GraphForBBfit->GetN(); ip ++) {
    protonV0GraphForBBfit->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonV0GraphForBBfit->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  for(Int_t ip = 0; ip < protonV0plusTOFGraph->GetN(); ip ++) {
    protonV0plusTOFGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonV0plusTOFGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  for(Int_t ip = 0; ip < protonV0plusTOFGraphForBBfit->GetN(); ip ++) {
    protonV0plusTOFGraphForBBfit->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonV0plusTOFGraphForBBfit->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  //
  //
  TGraphErrors * graphAll = new TGraphErrors();
  TList * listColl = new TList();
  if (!useV0s) {
    //listColl->Add(protonTpcGraph);
    //listColl->Add(kaonTpcGraph);  
    listColl->Add(pionTpcGraph);
    listColl->Add(electronGraph);
    listColl->Add(protonTofGraph); 
    listColl->Add(kaonTofGraph);
    listColl->Add(pionTofGraph);
  }
  else {
    listColl->Add(pionV0GraphForBBfit);
    //listColl->Add(electronV0GraphForBBfit);//nschmidt2016 for PbPb was commented for pp 
    listColl->Add(protonV0GraphForBBfit);
    
    //listColl->Add(pionV0plusTOFGraphForBBfit);
    listColl->Add(electronV0plusTOFGraphForBBfit);//nschmidt2016 for PbPb was uncommented for pp
    //listColl->Add(protonV0plusTOFGraphForBBfit);
  }
  MergeGraphErrors(graphAll, listColl);
  graphAll->SetNameTitle("beamDataPoints","beamDataPoints");
  
  // Create graph out of TPC only, V0, V0+TOF and TPC+TOF to determine correction functions
  TList * listCollPr = new TList();
  if (useV0s) {
    listCollPr->Add(protonTpcGraph);
    //listCollPr->Add(protonTofGraph);
    
    TGraphErrors * protonV0GraphCut = new TGraphErrors(*protonV0Graph);
    protonV0GraphCut->SetName("protonV0GraphCut");
    for(Int_t ip = 0; ip < protonV0GraphCut->GetN(); ip ++) {
      Double_t mom = protonV0GraphCut->GetX()[ip] * AliPID::ParticleMass(AliPID::kProton);
      //Bool_t removePoint = mom >= 0.6 || mom < 0.35; // Will take TPC protons instead (statistics) for very low momenta
      Bool_t removePoint = mom < (isPPb ? 0.6 : 0.35); // Will take TPC protons instead (statistics) for very low momenta;
                                                       // for pPb low statistics for the moment, therefore, go a bit further with TPC protons
      
      
      if (removePoint) {
        protonV0GraphCut->RemovePoint(ip);
        ip--;
        continue;
      }
    }
    /*
    TGraphErrors * protonV0plusTOFGraphCut = new TGraphErrors(*protonV0plusTOFGraph);
    protonV0plusTOFGraphCut->SetName("protonV0plusTOFGraphCut");
    for(Int_t ip = 0; ip < protonV0plusTOFGraphCut->GetN(); ip ++) {
      Double_t mom = protonV0plusTOFGraphCut->GetX()[ip] * AliPID::ParticleMass(AliPID::kProton);
      Bool_t removePoint = mom < 0.6;
      
      if (removePoint) {
        protonV0plusTOFGraphCut->RemovePoint(ip);
        ip--;
        continue;
      }
    }*/
  
    listCollPr->Add(protonV0GraphCut);
    //listCollPr->Add(protonV0plusTOFGraphCut);
  }
  else {
    listCollPr->Add(protonTpcGraph);
    listCollPr->Add(protonTofGraph);
  }
  TGraphErrors * graphPrCombined = new TGraphErrors();
  MergeGraphErrors(graphPrCombined, listCollPr);
  graphPrCombined->SetMarkerStyle(22);
  graphPrCombined->SetMarkerColor(4);
  graphPrCombined->SetNameTitle("Graph_Protons_Combined", "Graph_Protons_Combined");
  
  TList * listCollPi = new TList();
  if (useV0s) {
    //listCollPi->Add(pionTpcGraph);
    //listCollPi->Add(pionTofGraph);
    
    TGraphErrors * pionV0GraphCut = new TGraphErrors(*pionV0Graph);
    pionV0GraphCut->SetName("pionV0GraphCut");
    for(Int_t ip = 0; ip < pionV0GraphCut->GetN(); ip ++) {
      //Double_t mom = pionV0GraphCut->GetX()[ip] * AliPID::ParticleMass(AliPID::kPion);
      Bool_t removePoint = kFALSE;//mom >= 0.4;
                           //|| mom < 0.2;//TODO NOW NEEDED?
      if (removePoint) {
        pionV0GraphCut->RemovePoint(ip);
        ip--;
        continue;
      }
    }
    /*
    TGraphErrors * pionV0plusTOFGraphCut = new TGraphErrors(*pionV0plusTOFGraph);
    pionV0plusTOFGraphCut->SetName("pionV0plusTOFGraphCut");
    for(Int_t ip = 0; ip < pionV0plusTOFGraphCut->GetN(); ip ++) {
      Double_t mom = pionV0plusTOFGraphCut->GetX()[ip] * AliPID::ParticleMass(AliPID::kPion);
      Bool_t removePoint = mom < 0.4;
      
      if (removePoint) {
        pionV0plusTOFGraphCut->RemovePoint(ip);
        ip--;
        continue;
      }
    }*/
  
    listCollPi->Add(pionV0GraphCut);
    //listCollPi->Add(pionV0plusTOFGraphCut);
  }
  else {
    listCollPi->Add(pionTpcGraph);
    listCollPi->Add(pionTofGraph);
  }
  TGraphErrors * graphPiCombined = new TGraphErrors();
  MergeGraphErrors(graphPiCombined, listCollPi);
  graphPiCombined->SetMarkerStyle(22);
  graphPiCombined->SetMarkerColor(3);
  graphPiCombined->SetNameTitle("Graph_Pions_Combined", "Graph_Pions_Combined");
  
  TList * listCollKa = new TList();
  listCollKa->Add(kaonTpcGraph);
  listCollKa->Add(kaonTofGraph);
  TGraphErrors * graphKaCombined = new TGraphErrors();
  MergeGraphErrors(graphKaCombined, listCollKa);
  graphKaCombined->SetMarkerStyle(22);
  graphKaCombined->SetMarkerColor(5);
  graphKaCombined->SetNameTitle("Graph_Kaons_Combined", "Graph_Kaons_Combined");
  
  TList * listCollEl = new TList();
  if (useV0s) {
    //listCollEl->Add(electronGraph);
    
    TGraphErrors * electronV0GraphCut = new TGraphErrors(*electronV0Graph);
    electronV0GraphCut->SetName("electronV0GraphCut");
    /*
    for(Int_t ip = 0; ip < electronV0GraphCut->GetN(); ip ++) {
      Double_t mom = electronV0GraphCut->GetX()[ip] * AliPID::ParticleMass(AliPID::kElectron);
      Bool_t removePoint = mom > 0.3;
      if (removePoint) {
        electronV0GraphCut->RemovePoint(ip);
        ip--;
        continue;
      }
    }
    
    TGraphErrors * electronV0plusTOFGraphCut = new TGraphErrors(*electronV0plusTOFGraph);
    electronV0plusTOFGraphCut->SetName("electronV0plusTOFGraphCut");
    for(Int_t ip = 0; ip < electronV0plusTOFGraphCut->GetN(); ip ++) {
      Double_t mom = electronV0plusTOFGraphCut->GetX()[ip] * AliPID::ParticleMass(AliPID::kElectron);
      Bool_t removePoint = mom <= 0.4;
      
      if (removePoint) {
        electronV0plusTOFGraphCut->RemovePoint(ip);
        ip--;
        continue;
      }
    }
    */
    listCollEl->Add(electronV0GraphCut);
    //listCollEl->Add(electronV0plusTOFGraphCut);
  }
  else {
    listCollEl->Add(electronGraph);
  }
  TGraphErrors * graphElCombined = new TGraphErrors();
  MergeGraphErrors(graphElCombined, listCollEl);
  graphElCombined->SetMarkerStyle(22);
  graphElCombined->SetMarkerColor(2);
  graphElCombined->SetNameTitle("Graph_Electrons_Combined", "Graph_Electrons_Combined");
  
  
  
  //
  //  return array with all graphs
  //
  TObjArray * arrGraphs = new TObjArray();
  arrGraphs->Add(graphAll);
  //  
  arrGraphs->Add(electronGraph);
  arrGraphs->Add(electronV0Graph);
  arrGraphs->Add(electronV0plusTOFGraph);
  arrGraphs->Add(electronV0GraphForBBfit);
  arrGraphs->Add(pionTpcGraph);
  arrGraphs->Add(pionTofGraph);
  arrGraphs->Add(pionV0Graph);
  arrGraphs->Add(pionV0plusTOFGraph);
  arrGraphs->Add(pionV0GraphForBBfit);
  arrGraphs->Add(kaonTpcGraph);
  arrGraphs->Add(kaonTofGraph);
  arrGraphs->Add(protonTpcGraph);
  arrGraphs->Add(protonTofGraph);
  arrGraphs->Add(protonV0Graph);
  arrGraphs->Add(protonV0plusTOFGraph);
  arrGraphs->Add(protonV0GraphForBBfit);
  //
  arrGraphs->Add(graphElCombined);
  arrGraphs->Add(graphPiCombined);
  arrGraphs->Add(graphKaCombined);
  arrGraphs->Add(graphPrCombined);
  //
  //
  
  canvasQAtpc->SaveAs("splines_QA_ResidualGraphsTPC.root");
  canvasQAtof->SaveAs("splines_QA_ResidualGraphsTOF.root");
  canvasQAv0->SaveAs("splines_QA_ResidualGraphsV0.root");
  canvasQAv0plusTOF->SaveAs("splines_QA_ResidualGraphsV0plusTOF.root");
  
  
  canvasQAv0DeDxPurityEl->SaveAs("V0_dEdx_purity_el.root");
  canvasQAv0DeDxPurityPi->SaveAs("V0_dEdx_purity_Pi.root");
  canvasQAv0DeDxPurityPr->SaveAs("V0_dEdx_purity_Pr.root");

  return arrGraphs;

}


//________________________________________________________________________
TObjArray * AliTPCcalibResidualPID::GetResidualGraphsMC(THnSparseF * histPidQA, const Char_t * /*system*/) {
  //
  // Extracts the residual graphs from THnSparse created from MC (NOT DATA!)
  //
  
  const TString momTitle = "#it{p}_{TPC} (GeV/#it{c})";
  const TString dEdxTitle = "d#it{E}/d#it{x} (arb. unit)";
  
  Int_t cutForFitting = 10;
  Double_t heightFractionForFittingRange = 0.1; 
  
  //
  THnSparse * hist = histPidQA;
  //
  TCanvas * canvasQAmc = new TCanvas("canvasQAmcResGraph","Control canvas for residual graphs (MC)",100,10,1380,800);
  canvasQAmc->Divide(2,2);
  
  for (Int_t i = 1; i <= 4; i++) {
    canvasQAmc->GetPad(i)->SetGrid(1, 1);
    canvasQAmc->GetPad(i)->SetLogz();
    canvasQAmc->GetPad(i)->SetLogx();
  }
  //
  // 1. select and fit electrons
  //
  // In the following, axis 3 will also be limited in range, although nsigma itself is not used. This is to avoid multiple counting of entries
  hist->GetAxis(2)->SetRangeUser(0,0); // electrons
  hist->GetAxis(3)->SetRangeUser(0,0); // electrons (used for nsigma and for the eta correction of the signal)
  TH2D * histElectron = hist->Projection(1,0);
  histElectron->SetName("histElectron");
  histElectron->RebinX(2);
  histElectron->GetXaxis()->SetRangeUser(0.15,8.0);
  TObjArray arr;
  FitSlicesY(histElectron, heightFractionForFittingRange, cutForFitting, "QNRL", &arr);
  TH1D * electronPoints = (TH1D *) arr.At(1)->Clone();
  //
  TGraphErrors * electronGraph = new TGraphErrors(electronPoints);
  electronGraph->SetName("electronGraphMC");
  for(Int_t ip = 0; ip < electronGraph->GetN(); ip ++) {
    Bool_t removePoint = electronGraph->GetY()[ip] < 10 || electronGraph->GetEY()[ip]/electronGraph->GetY()[ip] > 0.03
                         || electronGraph->GetX()[ip] < 0.15;// || electronGraph->GetX()[ip] > 3.1;
    if (removePoint) {
      electronGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAmc->cd(1);
  histElectron->SetTitle("electrons");
  histElectron->GetYaxis()->SetRangeUser(50, 120);
  histElectron->GetXaxis()->SetTitle(momTitle.Data());
  histElectron->GetYaxis()->SetTitle(dEdxTitle.Data());
  histElectron->GetXaxis()->SetMoreLogLabels(kTRUE);
  histElectron->GetXaxis()->SetNoExponent(kTRUE);
  histElectron->Draw("colz");
  electronPoints->SetMarkerStyle(24);
  electronPoints->Draw("same");
  //
  // 2. protons
  //
  hist->GetAxis(2)->SetRangeUser(3,3); // protons
  hist->GetAxis(3)->SetRangeUser(3,3); // protons (used for nsigma and for the eta correction of the signal)
  //
  TH2D * histProton = hist->Projection(1,0);
  histProton->SetName("histProton");
  histProton->GetXaxis()->SetRangeUser(0.15,8);
  histProton->GetYaxis()->SetRangeUser(10, hist->GetAxis(1)->GetBinUpEdge(hist->GetAxis(1)->GetNbins()));
  FitSlicesY(histProton, heightFractionForFittingRange, cutForFitting, "QNRL", &arr);
  TH1D * protonPoints = (TH1D *) arr.At(1)->Clone();
  //
  TGraphErrors * protonGraph = new TGraphErrors(protonPoints);
  protonGraph->SetName("protonGraphMC");
  for(Int_t ip = 0; ip < protonGraph->GetN(); ip ++) {
    Bool_t removePoint = protonGraph->GetY()[ip] < 10 || protonGraph->GetEY()[ip]/protonGraph->GetY()[ip] > 0.02 || 
                         protonGraph->GetX()[ip] < 0.15 || protonGraph->GetX()[ip] > 6;
    if (removePoint) {
      protonGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAmc->cd(2);
  histProton->SetTitle("protons");
  histProton->GetXaxis()->SetTitle(momTitle.Data());
  histProton->GetYaxis()->SetTitle(dEdxTitle.Data());
  histProton->GetXaxis()->SetMoreLogLabels(kTRUE);
  histProton->GetXaxis()->SetNoExponent(kTRUE);
  histProton->Draw("colz");
  protonPoints->SetMarkerStyle(20);
  protonPoints->Draw("same");
  
  //
  // 3. pions
  //
  hist->GetAxis(2)->SetRangeUser(1,1); // pions
  hist->GetAxis(3)->SetRangeUser(1,1); // pions (used for nsigma and for the eta correction of the signal)
  //
  TH2D * histPion = hist->Projection(1,0);
  histPion->SetName("histPion");
  histPion->GetXaxis()->SetRangeUser(0.15,50);
  FitSlicesY(histPion, heightFractionForFittingRange, cutForFitting, "QNRL", &arr);
  TH1D * pionPoints =  (TH1D *) arr.At(1)->Clone();
  //
  TGraphErrors * pionGraph = new TGraphErrors(pionPoints);
  pionGraph->SetName("pionGraphMC");
  for(Int_t ip = 0; ip < pionGraph->GetN(); ip ++) {
    Bool_t removePoint = pionGraph->GetY()[ip] < 10 || pionGraph->GetEY()[ip]/pionGraph->GetY()[ip] > 0.005
                         || pionGraph->GetX()[ip] < 0.15 || pionGraph->GetX()[ip] > 30;
    if (removePoint) {
      pionGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAmc->cd(3);
  histPion->GetYaxis()->SetRangeUser(30, 90);
  histPion->GetXaxis()->SetTitle(momTitle.Data());
  histPion->GetYaxis()->SetTitle(dEdxTitle.Data());
  histPion->GetXaxis()->SetMoreLogLabels(kTRUE);
  histPion->GetXaxis()->SetNoExponent(kTRUE);
  histPion->Draw("colz");
  histPion->SetTitle("pions");
  pionPoints->SetMarkerStyle(20);
  pionPoints->Draw("same");
  
  //
  // 4. kaons
  //
  hist->GetAxis(2)->SetRangeUser(2,2); // kaons
  hist->GetAxis(3)->SetRangeUser(2,2); // kaons (used for nsigma and for the eta correction of the signal)
  //
  TH2D * histKaon = hist->Projection(1,0);
  histKaon->SetName("histKaon");
  histKaon->GetXaxis()->SetRangeUser(0.15,8);
  FitSlicesY(histKaon, heightFractionForFittingRange, cutForFitting, "QNRL", &arr);
  TH1D * kaonPoints = (TH1D*) arr.At(1)->Clone();
  //
  TGraphErrors * kaonGraph = new TGraphErrors(kaonPoints);
  kaonGraph->SetName("kaonGraphMC");
  for(Int_t ip = 0; ip < kaonGraph->GetN(); ip ++) {
    Bool_t removePoint = kaonGraph->GetY()[ip] < 10 || kaonGraph->GetEY()[ip]/kaonGraph->GetY()[ip] > 0.02
                         || kaonGraph->GetX()[ip] < 0.15;
    if (removePoint) {
      kaonGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  canvasQAmc->cd(4);
  histKaon->SetTitle("kaons");
  histKaon->GetYaxis()->SetRangeUser(30, 500);
  histKaon->GetXaxis()->SetTitle(momTitle.Data());
  histKaon->GetYaxis()->SetTitle(dEdxTitle.Data());
  histKaon->GetXaxis()->SetMoreLogLabels(kTRUE);
  histKaon->GetXaxis()->SetNoExponent(kTRUE);
  histKaon->Draw("colz");
  kaonPoints->SetMarkerStyle(20);
  kaonPoints->Draw("same");
  
  
  
  hist->GetAxis(2)->SetRange(0,-1); // RESET RANGES
  hist->GetAxis(3)->SetRange(0,-1); // RESET RANGES
  
  
  // Clone graphs to have graphs with the same names as for the data case -> same subsequent functions can be used to determine
  // the correction functions and, finally, the response functions.
  // In addition: Remove some points, which are expected to deviate from the Bethe-Bloch curve, from the original graphs, so that
  // these graphs can be used for the BB fit
  TGraphErrors * graphPrCombined = new TGraphErrors(*protonGraph);
  graphPrCombined->SetMarkerStyle(22);
  graphPrCombined->SetMarkerColor(4);
  graphPrCombined->SetNameTitle("Graph_Protons_Combined", "Graph_Protons_Combined");
  
  TGraphErrors * graphPiCombined = new TGraphErrors(*pionGraph);
  graphPiCombined->SetMarkerStyle(22);
  graphPiCombined->SetMarkerColor(3);
  graphPiCombined->SetNameTitle("Graph_Pions_Combined", "Graph_Pions_Combined");
  
  TGraphErrors * graphKaCombined = new TGraphErrors(*kaonGraph);
  graphKaCombined->SetMarkerStyle(22);
  graphKaCombined->SetMarkerColor(5);
  graphKaCombined->SetNameTitle("Graph_Kaons_Combined", "Graph_Kaons_Combined");
  
  TGraphErrors * graphElCombined = new TGraphErrors(*electronGraph);
  graphElCombined->SetMarkerStyle(22);
  graphElCombined->SetMarkerColor(2);
  graphElCombined->SetNameTitle("Graph_Electrons_Combined", "Graph_Electrons_Combined");
  
  for(Int_t ip = 0; ip < electronGraph->GetN(); ip ++) {
    Bool_t removePoint = electronGraph->GetX()[ip] < 2.5;// || electronGraph->GetX()[ip] > 5.0;
    //Bool_t removePoint = electronGraph->GetX()[ip] < 1.2 || electronGraph->GetX()[ip] > 5.0;
    if (removePoint) {
      electronGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  for(Int_t ip = 0; ip < protonGraph->GetN(); ip ++) {
    Bool_t removePoint = protonGraph->GetX()[ip] < 0.9
                         || protonGraph->GetX()[ip] > 2.5;//3.5; //TODO NEW in order not to get into conflict with pions/kaons
    if (removePoint) {
      protonGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  for(Int_t ip = 0; ip < pionGraph->GetN(); ip ++) {
    Bool_t removePoint = pionGraph->GetX()[ip] < 2.1;//TODO APR18 0.8;
    if (removePoint) {
      pionGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  
  for(Int_t ip = 0; ip < kaonGraph->GetN(); ip ++) {
    Bool_t removePoint = kaonGraph->GetX()[ip] < 1.25 || kaonGraph->GetX()[ip] > 4.0;// < 0.7;// || kaonGraph->GetX()[ip] > 4;
    if (removePoint) {
      kaonGraph->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  //
  // rescale graphs to betaGamma
  //
  kaonGraph->SetMarkerStyle(22);
  for(Int_t ip = 0; ip < kaonGraph->GetN(); ip ++) {
    kaonGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
    kaonGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
  }
  for(Int_t ip = 0; ip < graphKaCombined->GetN(); ip ++) {
    graphKaCombined->GetX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
    graphKaCombined->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kKaon);
  }
  //
  electronGraph->SetMarkerStyle(22);
  electronGraph->SetMarkerColor(2);
  for(Int_t ip = 0; ip < electronGraph->GetN(); ip ++) {
    electronGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    electronGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  for(Int_t ip = 0; ip < graphElCombined->GetN(); ip ++) {
    graphElCombined->GetX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
    graphElCombined->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kElectron);
  }
  //
  pionGraph->SetMarkerStyle(22);
  pionGraph->SetMarkerColor(3);
  for(Int_t ip = 0; ip < pionGraph->GetN(); ip ++) {
    pionGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    pionGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  for(Int_t ip = 0; ip < graphPiCombined->GetN(); ip ++) {
    graphPiCombined->GetX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
    graphPiCombined->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kPion);
  }
  //
  protonGraph->SetMarkerStyle(22);
  protonGraph->SetMarkerColor(4);
  for(Int_t ip = 0; ip < protonGraph->GetN(); ip ++) {
    protonGraph->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    protonGraph->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  for(Int_t ip = 0; ip < graphPrCombined->GetN(); ip ++) {
    graphPrCombined->GetX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
    graphPrCombined->GetEX()[ip] /= AliPID::ParticleMass(AliPID::kProton);
  }
  //
  //
  
  // Store graphs without further restriction
  TGraphErrors * graphPrAll = new TGraphErrors(*graphPrCombined);
  graphPrAll->SetNameTitle("Protons_MC_all", "Protons_MC_all");
  
  TGraphErrors * graphPiAll = new TGraphErrors(*graphPiCombined);
  graphPiAll->SetNameTitle("Pions_MC_all", "Pions_MC_all");
  
  TGraphErrors * graphKaAll = new TGraphErrors(*graphKaCombined);
  graphKaAll->SetNameTitle("Kaons_MC_all", "Kaons_MC_all");
  
  TGraphErrors * graphElAll = new TGraphErrors(*graphElCombined);
  graphElAll->SetNameTitle("Electrons_MC_all", "Electrons_MC_all");
  
  //
  //
  TGraphErrors * graphAll = new TGraphErrors();
  TList * listColl = new TList();
  listColl->Add(electronGraph);
  listColl->Add(protonGraph); 
  listColl->Add(kaonGraph);
  listColl->Add(pionGraph);
  MergeGraphErrors(graphAll, listColl);
  graphAll->SetNameTitle("beamDataPoints","beamDataPoints");
  
  //
  //  return array with all graphs
  //
  TObjArray * arrGraphs = new TObjArray();
  arrGraphs->Add(graphAll);
  //  
  arrGraphs->Add(electronGraph);
  arrGraphs->Add(pionGraph);
  arrGraphs->Add(kaonGraph);
  arrGraphs->Add(protonGraph);
  //
  arrGraphs->Add(graphElAll);
  arrGraphs->Add(graphPiAll);
  arrGraphs->Add(graphKaAll);
  arrGraphs->Add(graphPrAll);  
  //
  arrGraphs->Add(graphElCombined);
  arrGraphs->Add(graphPiCombined);
  arrGraphs->Add(graphKaCombined);
  arrGraphs->Add(graphPrCombined);
  //
  //
  
  canvasQAmc->SaveAs("splines_QA_ResidualGraphsTPC.root");

  return arrGraphs;

}


//________________________________________________________________________
TObjArray * AliTPCcalibResidualPID::GetResponseFunctions(TF1* parametrisation, TObjArray* inputGraphs, const Char_t * type, const Char_t * period, const Char_t * pass, const Char_t * system, const Char_t* dedxtype) {
  //
  //
  
  Bool_t isMC = kFALSE;
  if (strcmp(type, "MC") == 0) {
    isMC = kTRUE;
  }
  else if (strcmp(type, "DATA") == 0) {
    isMC = kFALSE;
  }
  else {
    Printf("ERROR - GetResponseFunctions: Unknown type \"%s\" - must be \"MC\" or \"DATA\"!", type);
    
    return 0x0;
  }
  
  Bool_t isPbPb = kFALSE;
  Bool_t isPPb = kFALSE;
  if (strcmp(system, "PBPB") == 0) {
    Printf("PbPb detected - Applying special handling!");
    isPbPb = kTRUE;
  }
  else if (strcmp(system, "PPB") == 0 || strcmp(system, "PBP") == 0) {
    Printf("p-Pb/Pb-p detected - Applying special handling!");
    isPPb = kTRUE;
  }
  
  TCanvas* canvasQA[4] = { 0x0, };
  for (Int_t i = 0; i < 4; i++) {
    canvasQA[i] = new TCanvas(Form("canvasQAres_%d", i), "Control canvas for residual polynomial",100,10,1380,800);
    canvasQA[i]->SetGrid(1, 1);
    canvasQA[i]->SetLogx();
  }
  /*
  TCanvas * canvasQA = new TCanvas("canvasQAres","Control canvas for residual polynomial",100,10,1380,800);
  canvasQA->Divide(2,2);
  
  for (Int_t i = 1; i <= 4; i++) {
    canvasQA->GetPad(i)->SetGrid(1, 1);
  }
  */
  //
  // smallest needed (100MeV proton) bg = 0.1/0.938(=AliPID::ParticleMass(AliPID::kProton)) = 0.1, 
  // largest needed 100 GeV electron, bg = 100/0.000511(=AliPID::ParticleMass(AliPID::kElectron)) = 200 000
  //
  Double_t from = 0.1;
  Double_t to = 200000;
  
  TF1* funcBB = new TF1(*parametrisation);
  funcBB->SetLineColor(2);
  funcBB->SetLineWidth(1);
  
  TString fitFuncString = "[0] + [1] / x + [2] / (x**2) + [3] / (x**3) + [4] / (x**4) + [5] / (x**5) + [6] * x";
  TString fitFuncString2 = "pol6";
  //
  // 1. extract proton corrections
  //
  TGraphErrors *  graphProtonTPC = (TGraphErrors *) inputGraphs->FindObject("Graph_Protons_Combined");
  TGraphErrors *  graphProtonTPCfinal = new TGraphErrors(*graphProtonTPC);
  graphProtonTPCfinal->SetName("graphProtonTPCfinal");
  for(Int_t ip = 0; ip < graphProtonTPC->GetN(); ip ++) {
    graphProtonTPC->GetY()[ip] -= funcBB->Eval(graphProtonTPC->GetX()[ip]);
    graphProtonTPC->GetY()[ip] /= funcBB->Eval(graphProtonTPC->GetX()[ip]);
    graphProtonTPC->GetEY()[ip] /= funcBB->Eval(graphProtonTPC->GetX()[ip]);;
  }
  canvasQA[0]->cd();
  //canvasQA->cd(1);
  //gPad->SetLogx();
  graphProtonTPC->GetXaxis()->SetTitle("#beta#gamma");
  graphProtonTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  graphProtonTPC->GetXaxis()->SetNoExponent(kTRUE);
  graphProtonTPC->GetYaxis()->SetLabelSize(0.04);
  graphProtonTPC->GetYaxis()->SetTitleOffset(0.85);
  graphProtonTPC->GetXaxis()->SetTitleOffset(1.0);
  graphProtonTPC->GetYaxis()->SetTitle("(data - fit) / fit");
  graphProtonTPC->Draw("ap");
  
  if (graphProtonTPC->GetN() <= 0)  {
    Printf("ERROR - GetResponseFunctions: Proton graph has no data points!");
    
    return 0x0;
  }
  graphProtonTPC->Sort(); // Sort points along x. Next, the very first point will be used to determine the starting point of the correction function
  TF1 * funcCorrProton = new TF1("funcCorrProton", fitFuncString2.Data(), TMath::Max(graphProtonTPC->GetX()[0], 0.15),1.0); // TODO BEN was 0.18 - 0.85 //nschmidt2016 was ...,1.0); //was then 2.0
  graphProtonTPC->Fit(funcCorrProton, "QREX0M");
  //
  // 2. extract kaon corrections
  //
  TGraphErrors *  graphKaonTPC = (TGraphErrors *) inputGraphs->FindObject("Graph_Kaons_Combined");
  TGraphErrors *  graphKaonTPCfinal = new TGraphErrors(*graphKaonTPC);
  graphKaonTPCfinal->SetName("graphKaonTPCfinal");
  for(Int_t ip = 0; ip < graphKaonTPC->GetN(); ip ++) {
    graphKaonTPC->GetY()[ip] -= funcBB->Eval(graphKaonTPC->GetX()[ip]);
    graphKaonTPC->GetY()[ip] /= funcBB->Eval(graphKaonTPC->GetX()[ip]);
    graphKaonTPC->GetEY()[ip] /= funcBB->Eval(graphKaonTPC->GetX()[ip]);;
  }
  canvasQA[1]->cd();
  //canvasQA->cd(2);
  //gPad->SetLogx();
  graphKaonTPC->GetXaxis()->SetTitle("#beta#gamma");
  graphKaonTPC->GetYaxis()->SetTitle("(data - fit) / fit");
  graphKaonTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  graphKaonTPC->GetXaxis()->SetNoExponent(kTRUE);
  graphKaonTPC->GetYaxis()->SetLabelSize(0.04);
  graphKaonTPC->GetYaxis()->SetTitleOffset(0.85);
  graphKaonTPC->GetXaxis()->SetTitleOffset(1.0);
  graphKaonTPC->Draw("ap");
  
  Bool_t kaonGraphHasDataPoints = kTRUE;
  if (graphKaonTPC->GetN() <= 0)  {
    Printf("WARNING - GetResponseFunctions: Kaon graph has no data points! Kaon splines will just be pure BB without low-p correction!");
    kaonGraphHasDataPoints = kFALSE;
  }
  graphKaonTPC->Sort(); // Sort points along x. Next, the very first point will be used to determine the starting point of the correction function

  TF1 * funcCorrKaon = 0x0;

  if (!kaonGraphHasDataPoints) {
    // Just take dummy correction function, which is constant at zero (= no correction)
    funcCorrKaon = new TF1("funcCorrKaon", "0", 0.1, 1);
  }
  else {
    if (isMC) {
      funcCorrKaon = new TF1("funcCorrKaon", fitFuncString.Data(), TMath::Max(graphKaonTPC->GetX()[0], 0.25), isPPb ? 3.44 : 1.981);
      graphKaonTPC->Fit(funcCorrKaon, "QREX0M");
    }
    else {
      // In case of data there are sometimes problems to fit the shape with one function (could describe the overall deviation,
      // but will not cover all the details).
      // Nevertheless, this shape (including most of the "details") can be fitted with the following approach with two functions
      TF1 * funcCorrKaon1 = new TF1("funcCorrKaon1", fitFuncString.Data(),
                                    TMath::Max(graphKaonTPC->GetX()[0], 0.25), 1.981); 
      graphKaonTPC->Fit(funcCorrKaon1, "QREX0M", "same", TMath::Max(graphKaonTPC->GetX()[0], 0.25), 1.0);

      TF1 * funcCorrKaon2 = new TF1("funcCorrKaon2", fitFuncString2.Data(), TMath::Max(graphKaonTPC->GetX()[0], 0.25),  1.981);
      graphKaonTPC->Fit(funcCorrKaon2, "QREX0M", "same", (isMC ? 1.981 : 1.0), 1.981);

      funcCorrKaon = new TF1("funcCorrKaon",
                             "funcCorrKaon1 * 0.5*(1.+TMath::Erf((1 - x) / 0.1)) + ([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x) * 0.5*(1.+TMath::Erf((x - 1) / 0.1))",
                             TMath::Max(graphKaonTPC->GetX()[0], 0.25), 1.981);

      for (Int_t i = funcCorrKaon1->GetNpar(), j = 0; i < funcCorrKaon1->GetNpar() + funcCorrKaon2->GetNpar(); i++, j++) {
        funcCorrKaon->SetParameter(j, funcCorrKaon2->GetParameter(j));
      }
      funcCorrKaon->SetLineColor(kRed);
      funcCorrKaon->GetHistogram()->DrawClone("csame");
      //funcCorrKaon->Draw("same");
    }
    /*TODO
    TF1 * funcCorrKaon = new TF1("funcCorrKaon", fitFuncString.Data(),//TODO BEN was fitFuncString2
                                 TMath::Max(graphKaonTPC->GetX()[0], 0.25), (isMC ? 1.981 : 1.0)); //TODO BEN was 0.79 for data and 1.45 for MC
    graphKaonTPC->Fit(funcCorrKaon, "QREX0M");
    */
  }
  //
  // 3. extract pion corrections
  //
  TGraphErrors *  graphPionTPC = (TGraphErrors *)  inputGraphs->FindObject("Graph_Pions_Combined");
  TGraphErrors *  graphPionTPCfinal = new TGraphErrors(*graphPionTPC);
  graphPionTPCfinal->SetName("graphPionTPCfinal");
  for(Int_t ip = 0; ip < graphPionTPC->GetN(); ip ++) {
    graphPionTPC->GetY()[ip] -= funcBB->Eval(graphPionTPC->GetX()[ip]);
    graphPionTPC->GetY()[ip] /= funcBB->Eval(graphPionTPC->GetX()[ip]);
    graphPionTPC->GetEY()[ip] /= funcBB->Eval(graphPionTPC->GetX()[ip]);
  }
  canvasQA[2]->cd();
  //canvasQA->cd(3);
  //gPad->SetLogx();
  graphPionTPC->GetXaxis()->SetTitle("#beta#gamma");
  graphPionTPC->GetYaxis()->SetTitle("(data - fit) / fit");
  graphPionTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  graphPionTPC->GetXaxis()->SetNoExponent(kTRUE);
  graphPionTPC->GetYaxis()->SetLabelSize(0.04);
  graphPionTPC->GetYaxis()->SetTitleOffset(0.85);
  graphPionTPC->GetXaxis()->SetTitleOffset(1.0);
  graphPionTPC->Draw("ap");
  
  if (graphPionTPC->GetN() <= 0)  {
    Printf("ERROR - GetResponseFunctions: Pion graph has no data points!");
    
    return 0x0;
  }
  graphPionTPC->Sort(); // Sort points along x. Next, the very first point will be used to determine the starting point of the correction function
  // In case of PbPb (data, not MC) only correct down to 0.2 GeV/c; otherwise: contamination due to electrons
  //~ TF1 * funcCorrPion = new TF1("funcCorrPion", fitFuncString.Data(), ((isPbPb && !isMC) ? (0.2 / AliPID::ParticleMass(AliPID::kPion)) : graphPionTPC->GetX()[0]), isMC ? (isPPb ? 12.8 : (isPbPb ? 12.8 : 7.1)) : (isPPb ? 7.68 : 7.1)); 
  TF1 * funcCorrPion = new TF1("funcCorrPion", fitFuncString.Data(), ((isPbPb && !isMC) ? (0.2 / AliPID::ParticleMass(AliPID::kPion)) : graphPionTPC->GetX()[0]), isMC ? (isPPb ? 12.8 : (isPbPb ? 12.8 : 7.1)) : (isPPb ? 7.68 : (isPbPb ? 15.5 : 7.1))); //nschmidt2016 was 12.8 : 7.1 // was then 27.5
  graphPionTPC->Fit(funcCorrPion, "QREX0M");
  //
  // 4. extract electron corrections
  //
  TGraphErrors *  graphElectronTPC = (TGraphErrors *)  inputGraphs->FindObject("Graph_Electrons_Combined");
  TGraphErrors *  graphElectronTPCfinal = new TGraphErrors(*graphElectronTPC);
  graphElectronTPCfinal->SetName("graphElectronTPCfinal");
  for(Int_t ip = 0; ip < graphElectronTPC->GetN(); ip ++) {
    graphElectronTPC->GetY()[ip] -= funcBB->Eval(graphElectronTPC->GetX()[ip]);
    graphElectronTPC->GetY()[ip] /= funcBB->Eval(graphElectronTPC->GetX()[ip]);
    graphElectronTPC->GetEY()[ip] /= funcBB->Eval(graphElectronTPC->GetX()[ip]);;
  }
  canvasQA[3]->cd();
  //canvasQA->cd(4);
  //gPad->SetLogx();
  graphElectronTPC->GetXaxis()->SetTitle("#beta#gamma");
  graphElectronTPC->GetYaxis()->SetTitle("(data - fit) / fit");
  graphElectronTPC->GetXaxis()->SetMoreLogLabels(kTRUE);
  graphElectronTPC->GetXaxis()->SetNoExponent(kTRUE);
  graphElectronTPC->GetYaxis()->SetLabelSize(0.04);
  graphElectronTPC->GetYaxis()->SetTitleOffset(0.85);
  graphElectronTPC->GetXaxis()->SetTitleOffset(1.0);
  graphElectronTPC->Draw("ap");
  
  if (graphElectronTPC->GetN() <= 0)  {
    Printf("ERROR - GetResponseFunctions: Electron graph has no data points!");
    
    return 0x0;
  }
  graphElectronTPC->Sort(); // Sort points along x. Next, the very first point will be used to determine the starting point of the correction function
  // In case of PbPb (data, not MC) only correct down to 0.2 GeV/c; otherwise: contamination due to pions
  TF1 * funcCorrElectron = new TF1("funcCorrElectron", fitFuncString.Data(), (!isMC && isPbPb ? (0.2 / AliPID::ParticleMass(AliPID::kElectron)) :graphElectronTPC->GetX()[0]), (isMC ? 3565 : (isPPb ? 2900 : (isPbPb ? 2500 : 1920/*970*/))));// TODO was 1800 for pp data
  // NOTE: For data, the results are almost the same for fitFuncString and fitFuncString2. Maybe, fitFuncString2 is slightly better.
  graphElectronTPC->Fit(funcCorrElectron, "QREX0M");
  //
  // EXTRACT GRAPHS AND PUT THEM TO THE TOBJARRAY
  //
  const Int_t nBins = 500;
  Double_t xBetaGamma[nBins];
  Double_t yProton[nBins];
  Double_t yKaon[nBins];
  Double_t yPion[nBins];
  Double_t yElectron[nBins];
  Double_t yDefault[nBins];
  //
  // 
  //
  xBetaGamma[0] = from;
  Double_t factor = pow(to/from, 1./nBins);
  
  //
  for(Int_t kk = 0; kk < nBins; kk++) {
    if (kk > 0) xBetaGamma[kk] = factor * xBetaGamma[kk-1];
    yProton[kk] =  funcBB->Eval(xBetaGamma[kk])/50.;
    yPion[kk] =  funcBB->Eval(xBetaGamma[kk])/50.;
    yKaon[kk] =  funcBB->Eval(xBetaGamma[kk])/50.;
    yElectron[kk] =  funcBB->Eval(xBetaGamma[kk])/50.;
    yDefault[kk] = funcBB->Eval(xBetaGamma[kk])/50.;

    // Added by Ben
    Double_t widthFactor = 0.020;
    Double_t smoothProton   = 0.5*(TMath::Erf((funcCorrProton->GetXmax()-xBetaGamma[kk])*AliPID::ParticleMass(AliPID::kProton)/widthFactor) + 1);
    Double_t smoothKaon     = 0.5*(TMath::Erf((funcCorrKaon->GetXmax()-xBetaGamma[kk])*AliPID::ParticleMass(AliPID::kKaon)/widthFactor) + 1);
    Double_t smoothPion     = 0.5*(TMath::Erf((funcCorrPion->GetXmax()-xBetaGamma[kk])*AliPID::ParticleMass(AliPID::kPion)/widthFactor) + 1);
    Double_t smoothElectron = 0.5*(TMath::Erf((funcCorrElectron->GetXmax()-xBetaGamma[kk])*AliPID::ParticleMass(AliPID::kElectron)/widthFactor) + 1);

    if (xBetaGamma[kk] > funcCorrProton->GetXmax())
      yProton[kk] *= (1 + smoothProton*funcCorrProton->Eval(funcCorrProton->GetXmax()));
    // Correction is so large that one runs into trouble at low bg for Protons.
    // The deviation is smaller if one takes the lower bound of the correction function
    else if (xBetaGamma[kk] < funcCorrProton->GetXmin())
      yProton[kk] *= (1 + smoothProton*funcCorrProton->Eval(funcCorrProton->GetXmin()));
    else
      yProton[kk] *= (1 + smoothProton*funcCorrProton->Eval(xBetaGamma[kk]));
    
    if (xBetaGamma[kk] > funcCorrKaon->GetXmax())
      yKaon[kk] *= (1 + smoothKaon*funcCorrKaon->Eval(funcCorrKaon->GetXmax()));
    // Correction is so large that one runs into trouble at low bg for Kaons.
    // The deviation is smaller if one takes the lower bound of the correction function
    else if (xBetaGamma[kk] < funcCorrKaon->GetXmin())
      yKaon[kk] *= (1 + smoothKaon*funcCorrKaon->Eval(funcCorrKaon->GetXmin()));
    else
      yKaon[kk] *= (1 + smoothKaon*funcCorrKaon->Eval(xBetaGamma[kk]));

    if (xBetaGamma[kk] > funcCorrElectron->GetXmax())
      yElectron[kk] *= (1 + smoothElectron*funcCorrElectron->Eval(funcCorrElectron->GetXmax()));
    else if (xBetaGamma[kk] < funcCorrElectron->GetXmin())
      yElectron[kk] *= (1 + smoothElectron*funcCorrElectron->Eval(funcCorrElectron->GetXmin()));
    else
      yElectron[kk] *= (1 + smoothElectron*funcCorrElectron->Eval(xBetaGamma[kk]));
    
    // Only true for LHC10d.pass?: Seems not to be needed because any (very small) deviations are most likely due to electron contamination(BEN)
    if (xBetaGamma[kk] > funcCorrPion->GetXmax())
      yPion[kk] *= (1 + smoothPion*funcCorrPion->Eval(funcCorrPion->GetXmax()));
    else if (xBetaGamma[kk] < funcCorrPion->GetXmin())
      yPion[kk] *= (1 + smoothPion*funcCorrPion->Eval(funcCorrPion->GetXmin()));
    else
      yPion[kk] *= (1 + smoothPion*funcCorrPion->Eval(xBetaGamma[kk]));

    
    /* Removed by Ben
    Double_t smoothProton = 0.5*(TMath::Erf((funcCorrProton->GetXmax()-xBetaGamma[kk])/0.002) + 1);
    Double_t smoothKaon   = 0.5*(TMath::Erf((funcCorrKaon->GetXmax()-xBetaGamma[kk])/0.002) + 1);
    Double_t smoothPion   = 0.5*(TMath::Erf((funcCorrPion->GetXmax()-xBetaGamma[kk])/0.002) + 1);

    if (xBetaGamma[kk] < funcCorrProton->GetXmax() && xBetaGamma[kk] > funcCorrProton->GetXmin()) yProton[kk] *= (1 + smoothProton*funcCorrProton->Eval(xBetaGamma[kk])); 
    if (xBetaGamma[kk] < funcCorrKaon->GetXmax() && xBetaGamma[kk] > funcCorrKaon->GetXmin()) yKaon[kk] *= (1 + smoothKaon*funcCorrKaon->Eval(xBetaGamma[kk])); 
    if (xBetaGamma[kk] < funcCorrPion->GetXmax() && xBetaGamma[kk] > funcCorrPion->GetXmin()) yPion[kk] *= (1 + smoothPion*funcCorrPion->Eval(xBetaGamma[kk])); 
    //if (xBetaGamma[kk] < funcCorrElectron->GetXmax()) yElectron[kk] *= (1 + funcCorrElectron->Eval(xBetaGamma[kk])); 
    //
    if (xBetaGamma[kk] < funcCorrProton->GetXmin()) yProton[kk] *= (1 + funcCorrProton->Eval(funcCorrProton->GetXmin())); 
    if (xBetaGamma[kk] < funcCorrKaon->GetXmin()) yKaon[kk] *= (1 + funcCorrKaon->Eval(funcCorrKaon->GetXmin())); 
    if (xBetaGamma[kk] < funcCorrPion->GetXmin()) yPion[kk] *= (1 + funcCorrPion->Eval(funcCorrPion->GetXmin())); 
    */
  }
  //
  TGraph * graphProton = new TGraph(nBins, xBetaGamma, yProton);
  TGraph * graphElectron = new TGraph(nBins, xBetaGamma, yElectron);
  TGraph * graphKaon = new TGraph(nBins, xBetaGamma, yKaon);
  TGraph * graphPion = new TGraph(nBins, xBetaGamma, yPion);
  TGraph * graphDefault = new TGraph(nBins, xBetaGamma, yDefault);
  //
  //
  TSpline3 * splineProton = new TSpline3("splineProton", graphProton);
  TSpline3 * splineElectron = new TSpline3("splineElectron", graphElectron);
  TSpline3 * splineKaon = new TSpline3("splineKaon", graphKaon);
  TSpline3 * splinePion = new TSpline3("splinePion", graphPion);
  TSpline3 * splineDefault = new TSpline3("splineDefault", graphDefault);
  //
  //
  splineProton->SetNameTitle(Form("TSPLINE3_%s_PROTON_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype),
                             Form("TSPLINE3_%s_PROTON_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype));
  splineElectron->SetNameTitle(Form("TSPLINE3_%s_ELECTRON_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype),
                               Form("TSPLINE3_%s_ELECTRON_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype));
  splineKaon->SetNameTitle(Form("TSPLINE3_%s_KAON_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype),
                           Form("TSPLINE3_%s_KAON_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype));
  splinePion->SetNameTitle(Form("TSPLINE3_%s_PION_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype),
                           Form("TSPLINE3_%s_PION_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype));
  splineDefault->SetNameTitle(Form("TSPLINE3_%s_ALL_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype),
                              Form("TSPLINE3_%s_ALL_%s_%s_%s_%sMEAN",type,period,pass,system,dedxtype));
  //
  TObjArray * arrResponse = new TObjArray();
  arrResponse->SetName("TPCpidResponseFunctions");
  arrResponse->AddLast(splineProton);
  arrResponse->AddLast(splineElectron);
  arrResponse->AddLast(splineKaon);
  arrResponse->AddLast(splinePion);
  arrResponse->AddLast(splineDefault);
  //
  
  // Draw deviation from final results
  for(Int_t ip = 0; ip < graphProtonTPCfinal->GetN(); ip ++) {
    graphProtonTPCfinal->GetY()[ip] -= 50. * splineProton->Eval(graphProtonTPCfinal->GetX()[ip]);
    graphProtonTPCfinal->GetY()[ip] /= 50. * splineProton->Eval(graphProtonTPCfinal->GetX()[ip]);
    graphProtonTPCfinal->GetEY()[ip] /= 50. * splineProton->Eval(graphProtonTPCfinal->GetX()[ip]);;
  }
  canvasQA[0]->cd();
  //canvasQA->cd(1);
  graphProtonTPCfinal->SetMarkerStyle(26);
  graphProtonTPCfinal->Draw("psame");
  
  for(Int_t ip = 0; ip < graphKaonTPCfinal->GetN(); ip ++) {
    graphKaonTPCfinal->GetY()[ip] -= 50. * splineKaon->Eval(graphKaonTPCfinal->GetX()[ip]);
    graphKaonTPCfinal->GetY()[ip] /= 50. * splineKaon->Eval(graphKaonTPCfinal->GetX()[ip]);
    graphKaonTPCfinal->GetEY()[ip] /= 50. * splineKaon->Eval(graphKaonTPCfinal->GetX()[ip]);;
  }
  canvasQA[1]->cd();
  //canvasQA->cd(2);
  graphKaonTPCfinal->SetMarkerStyle(26);
  graphKaonTPCfinal->Draw("psame");
  
  for(Int_t ip = 0; ip < graphPionTPCfinal->GetN(); ip ++) {
    graphPionTPCfinal->GetY()[ip] -= 50. * splinePion->Eval(graphPionTPCfinal->GetX()[ip]);
    graphPionTPCfinal->GetY()[ip] /= 50. * splinePion->Eval(graphPionTPCfinal->GetX()[ip]);
    graphPionTPCfinal->GetEY()[ip] /= 50. * splinePion->Eval(graphPionTPCfinal->GetX()[ip]);;
  }
  canvasQA[2]->cd();
  //canvasQA->cd(3);
  graphPionTPCfinal->SetMarkerStyle(26);
  graphPionTPCfinal->Draw("psame");
  
  for(Int_t ip = 0; ip < graphElectronTPCfinal->GetN(); ip ++) {
    graphElectronTPCfinal->GetY()[ip] -= 50. * splineElectron->Eval(graphElectronTPCfinal->GetX()[ip]);
    graphElectronTPCfinal->GetY()[ip] /= 50. * splineElectron->Eval(graphElectronTPCfinal->GetX()[ip]);
    graphElectronTPCfinal->GetEY()[ip] /= 50. * splineElectron->Eval(graphElectronTPCfinal->GetX()[ip]);;
  }
  canvasQA[3]->cd();
  //canvasQA->cd(4);
  graphElectronTPCfinal->SetMarkerStyle(26);
  graphElectronTPCfinal->Draw("psame");
  
  TFile* fSave = TFile::Open("splines_QA_ResidualPolynomials.root", "RECREATE");
  fSave->cd();
  for (Int_t i = 0; i < 4; i++){
    canvasQA[i]->Write();  
    canvasQA[i]->SaveAs(Form("polynomials%d.pdf",i));  
  }
  fSave->Close();
  //canvasQA->SaveAs("splines_QA_ResidualPolynomials.root");
  
  delete funcBB;
  return arrResponse;
}

TF1* AliTPCcalibResidualPID::SetUpFitFunction(const Double_t* initialParameters, AliTPCcalibResidualPID::FitType fitType, Float_t from, Float_t to, Bool_t isPPb, Bool_t isMC, Double_t* parametersBBForward) {
  
  const Int_t nPar = 6;
  
  parametersBBForward = new Double_t[nPar];
  
  TF1* funcBB;

  if (fitType == AliTPCcalibResidualPID::kSaturatedLund) {
    printf("Fit function: Saturated Lund\n");
    
    funcBB = new TF1("SaturatedLund", SaturatedLund, from, to, nPar);
    //Double_t parametersBB[nPar] = {34.0446, 8.42221, 4.16724, 1.29473, 80.6663, 0}; //Xianguos values
    //Double_t parametersBB[nPar] = {35.5, 8.7, 2.0, 1.09, 75.6, 0}; // No saturation
    Double_t parametersBB[nPar] = {61.0, 8.7, 0.1, 0.85, 113.4, -38}; // Yields reasonable results for data and MC ~ all periods
    
    if (isPPb && !isMC) {
      parametersBB[0] = 51.6;
      parametersBB[1] = 9.7;
      parametersBB[2] = 1.62;
      parametersBB[3] = 0.99;
      parametersBB[4] = 104.4;
      parametersBB[5] = -27.0;
    }
    if (isMC) {
      parametersBB[0] = 41.4;
      parametersBB[1] = 8.6;
      parametersBB[2] = 2.2;
      parametersBB[3] = 0.92;
      parametersBB[4] = 90.4;
      parametersBB[5] = -20.0;
    }
    funcBB->SetParameters(parametersBB);
    
    Double_t parameterErrorsBB[nPar] = {5., 0.5, 0.2, 0.05, 10, 10};
    funcBB->SetParErrors(parameterErrorsBB);
    
    for (Int_t i = 0; i < nPar; i++)
      parametersBBForward[i] = parametersBB[i];
  }
  else if (fitType == AliTPCcalibResidualPID::kLund) {
    printf("Fit function: Lund\n");
    printf("****WARNING: Fit settings not tuned for this fit function, only for saturated Lund!\n");

    funcBB = new TF1("SaturatedLund", SaturatedLund, from, to, nPar);
    //Double_t parametersBB[nPar] = {34.0446, 8.42221, 4.16724, 1.29473, 80.6663, 0}; //Xianguos values
    //Double_t parametersBB[nPar] = {35.5, 8.7, 2.0, 1.09, 75.6, 0}; // No saturation
    Double_t parametersBB[nPar] = {35.0, 7., 1.86, 1., 75., 0}; // Yields reasonable results for data and MC ~ all periods
    
    funcBB->SetParameters(parametersBB);
    
    Double_t parameterErrorsBB[nPar] = {5., 0.5, 0.2, 0.05, 10, 10};
    funcBB->SetParErrors(parameterErrorsBB);
    
    funcBB->FixParameter(5, 0.0); // No saturation
    
    for (Int_t i = 0; i < nPar; i++)
      parametersBBForward[i] = parametersBB[i];
  }
  else if (fitType == kAlephWithAdditionalParam) {
    printf("Fit function: Aleph with additional parameter\n");
    printf("****WARNING: Fit settings not tuned for this fit function, only for saturated Lund!\n");
    
    // NOTE: This form is equivalent to the original form, but with parameter [2] redefined for better numerical stability.
    // The additional parameter [5] has been introduced later and is unity originally. It seems not to be needed and is, thus,
    // fixed to unity
    funcBB = new TF1("funcModifiedAleph", Aleph, from, to, nPar);
    
    Double_t parametersBB[nPar] = {1.6,20.3,-10,2.6,2.3, 0.02};
    //OLD with different sign for [3]: Double_t parametersBB[nPar] = {0.0762*50.3,10.632,TMath::Log(1.34e-05),1.863,1.948, 1};
    funcBB->SetParameters(parametersBB);
    
    for (Int_t i = 0; i < nPar; i++)
      parametersBBForward[i] = parametersBB[i];
  }
  else if (fitType == AliTPCcalibResidualPID::kAleph) {
    printf("Fit function: Aleph\n");
    printf("****WARNING: Fit settings not tuned for this fit function, only for saturated Lund!\n");
    
    // NOTE: This form is equivalent to the original form, but with parameter [2] redefined for better numerical stability.
    // The additional parameter [5] has been introduced later and is unity originally. It seems not to be needed and is, thus,
    // fixed to unity
    funcBB = new TF1("funcAleph",Aleph, from, to, nPar);
                     //OLD"[0]*([1]*TMath::Power(TMath::Sqrt(1 + x*x)/x , [3]) - [5] - TMath::Power(TMath::Sqrt(1 + x*x)/x , [3])*TMath::Log(TMath::Exp([2]) + 1/TMath::Power(x, [4])))", from, to); 
    //TEST Double_t parametersBB[nPar] = {1.2, 26.0, -30.0, -2.15, 5.6, 1};
    Double_t parametersBB[nPar] = {1.25, 27.5, -29.0, 2.2, 5.2, 1};
    
    // For [5] floating: Double_t parametersBB[nPar] = {2.42,15.2,-16,-2.24,2.8, 0.057};
    //OLD with different sign for [3] and [5] not fix: Double_t parametersBB[nPar] = {2.6,14.3,-15,2.2,2.7, 0.06};
    //OLD with different sign for [3]: Double_t parametersBB[nPar] = {0.0762*50.3,10.632,TMath::Log(1.34e-05),1.863,1.948, 1};
    funcBB->SetParameters(parametersBB);
    
    // Fix parameter 5 to original value of unity
    funcBB->FixParameter(5, 1); 
    
    for (Int_t i = 0; i < nPar; i++)
      parametersBBForward[i] = parametersBB[i];
  }
  else if (fitType == AliTPCcalibResidualPID::kAlephExternal) {
    printf("Fit function: Aleph from AliExternalTrackParam\n");
    printf("****WARNING: Fit settings not tuned for this fit function, only for saturated Lund!\n");

    // NOTE: This form is equivalent to the original form, but with parameter [2] redefined for better numerical stability.
    // The additional parameter [5] has been introduced later and is unity originally. It seems not to be needed and is, thus,
    // fixed to unity
    funcBB = new TF1("funcAleph",
                     "AliExternalTrackParam::BetheBlochAleph(x,[0],[1],[2],[3],[4])", from, to);
    Double_t parametersBB[nPar] = {0.0851148*50, 9.25771, 2.6558e-05, 2.32742, 1.83039, 0};

    funcBB->SetParameters(parametersBB);

    for (Int_t i = 0; i < nPar; i++)
      parametersBBForward[i] = parametersBB[i];
  }
  else {
    printf("Error: fitType not supported!\n");
    return 0x0;
  }
  
  funcBB->SetLineColor(2);
  funcBB->SetLineWidth(1);
  
  // Override initial parameters, if user provides some
  if (initialParameters) {
    printf("Setting user initial parameters!\n");
    funcBB->SetParameters(initialParameters);
  }
  
  return funcBB;
}

void AliTPCcalibResidualPID::SetUpInputGraph(TGraphErrors* graphAll, Bool_t isMC, Bool_t useV0s) {

  if (isMC) {
    for(Int_t ip = 0; ip < graphAll->GetN(); ip++) {
      graphAll->SetPointError(ip, 0, 0);
    }
  }
  else if (useV0s) {
    // Increase weight in plateau by reducing errors. Neccessary because errors are large there compared to other data points
    for(Int_t ip = 0; ip < graphAll->GetN(); ip++) {
      if (graphAll->GetX()[ip] >= 2500) { 
        graphAll->SetPointError(ip, 0, graphAll->GetEY()[ip] / 6.0); 
      }
      // Same for pions in the very rel. rise
      if (graphAll->GetX()[ip] >= 25 && graphAll->GetX()[ip] < 1000) {
        graphAll->SetPointError(ip, 0, graphAll->GetEY()[ip] / 6.0); 
      }
      // Same for protons in the MIP region
      if (graphAll->GetX()[ip] >= 2 && graphAll->GetX()[ip] < 4) {
        graphAll->SetPointError(ip, 0, graphAll->GetEY()[ip] / 6.0); 
      }
    }
  }
}

TCanvas* AliTPCcalibResidualPID::CreateBBCanvas(TObjArray* inputGraphs, Bool_t isMC, TF1* func) {
  
  TGraphErrors * graphAll = (TGraphErrors *) inputGraphs->FindObject("beamDataPoints");
  TCanvas * canvDelta_1 = new TCanvas("canvDelta_1","control histogram for Bethe-Bloch fit 1",100,10,1380,800);
  canvDelta_1->SetGrid(1, 1);
  canvDelta_1->SetLogx();  
  canvDelta_1->cd();
  TH1F *hBBdummy=new TH1F("hBBdummy","BB fit;#beta#gamma;#LTdE/dx#GT (arb. unit)",100,.8,1e4);
  hBBdummy->SetMinimum(45);
  hBBdummy->SetMaximum(120);
  hBBdummy->GetXaxis()->SetTitleOffset(1.1);
  hBBdummy->SetStats(kFALSE);
  hBBdummy->Draw();  
  
  graphAll->SetTitle("BB multi-graph fit");
  graphAll->SetMarkerStyle(22);
  graphAll->SetMarkerColor(kMagenta);
  graphAll->Draw("p");

  TLegend *leg=new TLegend(.2,0.3,0.36,isMC?.5:0.7);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);

  //overlay individual graphs
  TGraph *gr=0x0;
  if (isMC) {
    gr=(TGraph*)inputGraphs->FindObject("Protons_MC_all");
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"p (MC)","p");
    gr=(TGraph*)inputGraphs->FindObject("Pions_MC_all");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"#pi (MC)","p");
    gr=(TGraph*)inputGraphs->FindObject("Kaons_MC_all");
    gr->SetMarkerColor(kGreen);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"K (MC)","p");
    gr=(TGraph*)inputGraphs->FindObject("Electrons_MC_all");
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"e (MC)","p");
  }
  else {
    gr=(TGraph*)inputGraphs->FindObject("protonTpcGraph");
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(26);
    gr->Draw("p");
    leg->AddEntry(gr,"p (TPC)","p");
    gr=(TGraph*)inputGraphs->FindObject("protonTofGraph");
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(25);
    gr->Draw("p");
    leg->AddEntry(gr,"p (TOF)","p");
    gr=(TGraph*)inputGraphs->FindObject("protonV0Graph");//ForBBfit");
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"p (V0)","p");
    gr=(TGraph*)inputGraphs->FindObject("protonV0plusTOFGraph");//ForBBfit");
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(30);
    gr->Draw("p");
    leg->AddEntry(gr,"p (V0+TOF)","p");
    gr=(TGraph*)inputGraphs->FindObject("pionTpcGraph");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(26);
    gr->Draw("p");
    leg->AddEntry(gr,"#pi (TPC)","p");
    gr=(TGraph*)inputGraphs->FindObject("pionTofGraph");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(25);
    gr->Draw("p");
    leg->AddEntry(gr,"#pi (TOF)","p");
    gr=(TGraph*)inputGraphs->FindObject("pionV0Graph");//ForBBfit");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"#pi (V0)","p");
    gr=(TGraph*)inputGraphs->FindObject("pionV0plusTOFGraph");//ForBBfit");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(30);
    gr->Draw("p");
    leg->AddEntry(gr,"#pi (V0+TOF)","p");
    gr=(TGraph*)inputGraphs->FindObject("kaonTpcGraph");
    gr->SetMarkerColor(kGreen);
    gr->SetMarkerStyle(26);
    gr->Draw("p");
    leg->AddEntry(gr,"K (TPC)","p");
    gr=(TGraph*)inputGraphs->FindObject("kaonTofGraph");
    gr->SetMarkerColor(kGreen);
    gr->SetMarkerStyle(25);
    gr->Draw("p");
    leg->AddEntry(gr,"K (TOF)","p");
    gr=(TGraph*)inputGraphs->FindObject("electronGraph");
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(25);
    gr->Draw("p");
    leg->AddEntry(gr,"e","p");
    gr=(TGraph*)inputGraphs->FindObject("electronV0Graph");//ForBBfit");
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(24);
    gr->Draw("p");
    leg->AddEntry(gr,"e (V0)","p");
    gr=(TGraph*)inputGraphs->FindObject("electronV0plusTOFGraph");//ForBBfit");
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(30);
    gr->Draw("p");
    leg->AddEntry(gr,"e (V0+TOF)","p");
  }
  
  graphAll->GetXaxis()->SetTitle("#beta#gamma");
  graphAll->GetYaxis()->SetTitle("<dE/dx> (arb. unit)");
  graphAll->Draw("p");
  func->GetHistogram()->DrawClone("csame");
 
  leg->AddEntry(graphAll,"used","p");
  leg->AddEntry(func,"fit","l");
  leg->Draw("same"); 
  return canvDelta_1;
}

TCanvas* AliTPCcalibResidualPID::CreateResidualCanvas(TGraphErrors * graphAll, TF1* func) {
  TGraphErrors * graphDelta = new TGraphErrors(*graphAll);
  graphDelta->SetName("graphDelta");

  for(Int_t ip = 0; ip < graphDelta->GetN(); ip++) {
    graphDelta->GetY()[ip] -= func->Eval(graphDelta->GetX()[ip]);
    graphDelta->GetY()[ip] /= func->Eval(graphDelta->GetX()[ip]);
    graphDelta->GetEY()[ip] /= func->Eval(graphDelta->GetX()[ip]);
  }  
  
  TCanvas * canvDelta_2 = new TCanvas("canvDelta_2","control histogram for Bethe-Bloch fit 2",100,10,1380,800);
  canvDelta_2->SetGrid(1, 1);
  canvDelta_2->SetLogx();
  canvDelta_2->cd();   
  
  TH1F *hBBresdummy=new TH1F("hBBresdummy","residuals of BB fit;#beta#gamma;(data-fit)/fit",100,.8,1e4);
  hBBresdummy->SetMinimum(-0.04);
  hBBresdummy->SetMaximum(0.04);
  hBBresdummy->GetXaxis()->SetTitleOffset(1.1);
  hBBresdummy->GetYaxis()->SetTitleOffset(0.8);
  hBBresdummy->SetStats(kFALSE);
  hBBresdummy->Draw();

  graphDelta->SetTitle("residuals of BB fit to multi-graph");
  graphDelta->GetXaxis()->SetTitle("#beta#gamma");
  graphDelta->GetYaxis()->SetTitle("(data - fit) / fit");
  graphDelta->SetMarkerStyle(22);
  graphDelta->SetMarkerColor(4);
  graphDelta->Draw("p"); 
  
  return canvDelta_2;
}


//________________________________________________________________________
TF1* AliTPCcalibResidualPID::FitBB(TObjArray* inputGraphs, Bool_t isMC, Bool_t isPPb, const Bool_t useV0s, 
                                   const Double_t * initialParameters, AliTPCcalibResidualPID::FitType fitType) {
  //
  // Fit Bethe-Bloch parametrisation to data points (inputGraphs)
  //
  TGraphErrors * graphAll = (TGraphErrors *) inputGraphs->FindObject("beamDataPoints");
  //
  Float_t from = 0.9; //TODO ADJUST -> Very important
  Float_t to = graphAll->GetXaxis()->GetXmax() * 2;
  
  Double_t* parametersBBForward = 0x0;
  
  TVirtualFitter::SetMaxIterations(5e6);
  
  TF1* funcBB = SetUpFitFunction(initialParameters, fitType, from, to, isPPb, isMC, parametersBBForward);
  
  // In MC case: Remove errors from fit -> Some data points with extremely small errors otherwise completely dominate the fit,
  // but these points could still be wrong due to systematics (low p effects, e.g.)
  SetUpInputGraph(graphAll, isMC, useV0s);

  graphAll->Fit(funcBB, "REX0M");  
  funcBB->SetRange(from, to);
//   funcBB->GetParameters(parametersBBForward);
  
  TCanvas * canvDelta_1 = CreateBBCanvas(inputGraphs, isMC, funcBB);
  TCanvas * canvDelta_2 = CreateResidualCanvas(graphAll, funcBB);
   
  TString fitTypeName = (GetStringFitType(fitType)).ReplaceAll(" ","");
  TFile* fSave = TFile::Open(TString::Format("splines_QA_BetheBlochFit_%s.root",fitTypeName.Data()).Data(), "RECREATE");
  fSave->cd();
  canvDelta_1->Write();
  printf("Save canvas");
  canvDelta_1->SaveAs(TString::Format("bethebloch_%s.pdf",fitTypeName.Data()).Data());
  canvDelta_2->Write();
  canvDelta_2->SaveAs(TString::Format("betheblochresidual_%s.pdf",fitTypeName.Data()).Data());
  
  TString fitResults = "Fit results:\n";
  for (Int_t i = 0; i < funcBB->GetNpar(); i++) {
    fitResults.Append(Form("par%d:\t%f +- %f\n", i, funcBB->GetParameters()[i], funcBB->GetParErrors()[i]));
  }
  
  TNamed* settings = new TNamed(fitResults.Data(), fitResults.Data());
  settings->Write();
  fSave->Close();
  
  
  printf("\n\n%s", fitResults.Data());
  
  return funcBB;
}


//________________________________________________________________________
Double_t AliTPCcalibResidualPID::Lund(Double_t* xx, Double_t* par)
{
  // bg is beta*gamma
  const Double_t bg = xx[0];

  const Double_t beta2 = bg*bg / (1.0 + bg*bg);
  
  const Double_t a = par[0];
  const Double_t b = par[1];
  const Double_t c = par[2];
  const Double_t e = par[3];
  const Double_t f = par[4];
  
  const Double_t d = TMath::Exp(c*(a - f) / b);
 
  const Double_t powbg = TMath::Power(1.0 + bg, c);
 
  const Double_t value = a / TMath::Power(beta2,e) +
    b/c * TMath::Log(powbg / (1.0 + d*powbg));
    
  return value;
}


//________________________________________________________________________
Double_t AliTPCcalibResidualPID::SaturatedLund(Double_t* xx, Double_t* par)
{
  const Double_t qq = Lund(xx, par);
  return qq * TMath::Exp(par[5] / qq);
}

Double_t AliTPCcalibResidualPID::Aleph(Double_t* xx, Double_t* par) {
  const Double_t bg = xx[0];
  
  const Double_t a0 = par[0];
  const Double_t a1 = par[1];
  const Double_t a2 = par[2];
  const Double_t a3 = par[3];
  const Double_t a4 = par[4];
  const Double_t a5 = par[5];
  
  const Double_t beta = TMath::Sqrt(bg*bg / (1.0 + bg*bg));
  const Double_t powbetaa3 = TMath::Power(beta,a3);
  
  const Double_t value = a0/powbetaa3 * (a1 - a2 - a5 * powbetaa3 - TMath::Log(1.0 + TMath::Power(bg, -a4)*TMath::Exp(-a2)));

  return value;
}


//________________________________________________________________________
void AliTPCcalibResidualPID::FitSlicesY(TH2 *hist, Double_t heightFractionForRange, Int_t cutThreshold, TString fitOption, TObjArray *arr)
{
  //heightPercentageForRange
  // custom slices fit
  //

  if (!hist) return;
  if (!arr) return;

  // If old style is to be used
  /*
  hist->FitSlicesY(0, 0, -1, cutThreshold, fitOption.Data(), &arr);
  return;
  */
  
  
  arr->Clear();
  arr->SetOwner();
  arr->Expand(4);
  
  TAxis *axis=hist->GetXaxis();
  const TArrayD *bins = axis->GetXbins();
  TH1D** hList = new TH1D*[4];
  
  for (Int_t i = 0; i < 4; i++) {
    delete gDirectory->FindObject(Form("%s_%d", hist->GetName(), i));
    
    if (bins->fN == 0) {
      hList[i] = new TH1D(Form("%s_%d", hist->GetName(), i), i < 3 ? Form("Fitted value of par[%d]", i) : "Chi2/NDF",
                          hist->GetNbinsX(), axis->GetXmin(), axis->GetXmax());
    } else {
      hList[i] = new TH1D(Form("%s_%d", hist->GetName(), i), i < 3 ? Form("Fitted value of par[%d]", i) : "Chi2/NDF", hist->GetNbinsX(), bins->fArray);
    }
    
    (*arr)[i] = hList[i];
  }
  
  for (Int_t ibin=axis->GetFirst(); ibin<=axis->GetLast(); ++ibin){
    TH1 *h=hist->ProjectionY("_temp",ibin,ibin);
    if (!h)
      continue;
    
    if (h->GetEntries() < cutThreshold) {
      delete h;
      continue;
    }
    
    // Average around maximum bin -> Might catch some outliers
    Int_t maxBin = h->GetMaximumBin();
    Double_t maxVal = h->GetBinContent(maxBin);
    
    if (maxVal < 2) { // It could happen that all entries are in overflow/underflow; don't fit in this case
      delete h;
      continue;
    }
    
    UChar_t usedBins = 1;
    if (maxBin > 1) {
      maxVal += h->GetBinContent(maxBin - 1);
      usedBins++;
    }
    if (maxBin < h->GetNbinsX()) {
      maxVal += h->GetBinContent(maxBin + 1);
      usedBins++;
    }
    maxVal /= usedBins;
    
    Double_t thresholdFraction = heightFractionForRange * maxVal; 
    Int_t lowThrBin = TMath::Max(1, h->FindFirstBinAbove(thresholdFraction));
    Int_t highThrBin = TMath::Min(h->GetNbinsX(), h->FindLastBinAbove(thresholdFraction));
    
    Double_t lowThreshold = h->GetBinCenter(lowThrBin);
    Double_t highThreshold = h->GetBinCenter(highThrBin);

    TFitResultPtr res = h->Fit("gaus", Form("%sS", fitOption.Data()), "", lowThreshold, highThreshold);
    
    if ((Int_t)res == 0) {
      Int_t resBin = ibin;
      hList[0]->SetBinContent(resBin,res->GetParams()[0]);
      hList[0]->SetBinError(resBin,res->GetErrors()[0]);
      hList[1]->SetBinContent(resBin,res->GetParams()[1]);
      hList[1]->SetBinError(resBin,res->GetErrors()[1]);
      hList[2]->SetBinContent(resBin,res->GetParams()[2]);
      hList[2]->SetBinError(resBin,res->GetErrors()[2]);
      hList[3]->SetBinContent(resBin,res->Ndf()>0?res->Chi2()/res->Ndf():0);
    }
    
    delete h;
  }
}


//______________________________________________________________________________
Bool_t AliTPCcalibResidualPID::GetVertexIsOk(AliVEvent* event) const
{
  AliAODEvent* aod = 0x0;
  AliESDEvent* esd = 0x0;
  
  aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    esd = dynamic_cast<AliESDEvent*>(event);
    
    if (!esd) {
      AliError("Event seems to be neither AOD nor ESD!");
      return kFALSE;
    }
  }
    
  if (fIsPbpOrpPb) {
    const AliVVertex* trkVtx = (aod ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex()) :
                                      dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertex()));
      
    if (!trkVtx || trkVtx->GetNContributors() <= 0)
      return kFALSE;
      
    TString vtxTtl = trkVtx->GetTitle();
    if (!vtxTtl.Contains("VertexerTracks"))
      return kFALSE;
      
    Float_t zvtx = trkVtx->GetZ();
    const AliVVertex* spdVtx = (aod ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertexSPD()) :
                                      dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexSPD()));
    if (spdVtx->GetNContributors() <= 0)
      return kFALSE;
      
    Double_t cov[6] = {0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (spdVtx->IsFromVertexerZ() && (zRes > 0.25))
      return kFALSE;
      
    if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ()) > 0.5)
      return kFALSE;

    if (TMath::Abs(zvtx) > fZvtxCutEvent) //Default: 10 cm
      return kFALSE;
      
    return kTRUE;
  }
    
  
  // pp and PbPb
  const AliVVertex* primaryVertex = (aod ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex()) :
                                           dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexTracks()));
    
  if (!primaryVertex || primaryVertex->GetNContributors() <= 0)
    return kFALSE;
      
  if (TMath::Abs(primaryVertex->GetZ()) >= fZvtxCutEvent) //Default: 10 cm
    return kFALSE;
  
  return kTRUE;
}


//______________________________________________________________________________
void AliTPCcalibResidualPID::FillV0PIDlist(AliESDEvent* event)
{
  //
  // Fill the PID tag list
  //
  
  // If no event forwarded as parameter (default), cast current input event.
  // Dynamic cast to ESD events (DO NOTHING for AOD events)
  if (!event)
    event = dynamic_cast<AliESDEvent *>(InputEvent());
  
  // If this fails, do nothing
  if (!event)
    return;
  
  if (!fV0KineCuts) {
    AliError("V0KineCuts not available!");
    return;
  }
  
  TString beamType(event->GetBeamType());
  
  if (beamType.CompareTo("Pb-Pb") == 0 || beamType.CompareTo("A-A") == 0) {
    fV0KineCuts->SetMode(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb);
  }
  else {
    fV0KineCuts->SetMode(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPP); 
  }

  // V0 selection
  // set event
  fV0KineCuts->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Char_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++)
    fV0tags[i] = 0;
  
  fV0motherIndex = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++)
    fV0motherIndex[i] = -1;
  
  fV0motherPDG = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++)
    fV0motherPDG[i] = 0;
  
  fNumTagsStored = numTracks;
  
  
  
  
  
  // loop over V0 particles
  for (Int_t iv0 = 0; iv0 < event->GetNumberOfV0s(); iv0++) {
    AliESDv0* v0 = (AliESDv0*)event->GetV0(iv0);
 
    if (!v0)
      continue;
    
    // Reject onFly V0's <-> Only take V0's from offline V0 finder
    if (v0->GetOnFlyStatus())
      continue; 
  
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0 = 0, pdgP = 0, pdgN = 0;
    
    foundV0 = fV0KineCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if (!foundV0)
      continue;
    
    Int_t iTrackP = v0->GetPindex();  // positive track
    Int_t iTrackN = v0->GetNindex();  // negative track

    /*
    AliESDtrack* trackP = event->GetTrack(iTrackP);
    AliESDtrack* trackN = event->GetTrack(iTrackN);
    
    if (!trackP || !trackN)
      continue;
    
    
    
    Float_t xy = 999, z = 999;
    trackP->GetImpactParameters(xy, z);
    
    Bool_t reject = kFALSE;
    if (TMath::Abs(z) > 1)
      continue;
    
    trackN->GetImpactParameters(xy, z);
    if (TMath::Abs(z) > 1)
      continue;
    */
    
    
    // Fill the Object arrays
    // positive particles
    if (pdgP == -11) {
      fV0tags[iTrackP] = 14;
    }
    else if (pdgP == 211) {
      fV0tags[iTrackP] = 15;
    }
    else if(pdgP == 2212) {
      fV0tags[iTrackP] = 16;
    }
    
    fV0motherIndex[iTrackP] = iv0;
    fV0motherPDG[iTrackP] = pdgV0;

    // negative particles
    if( pdgN == 11){
      fV0tags[iTrackN] = -14;
    }
    else if( pdgN == -211){
      fV0tags[iTrackN] = -15;
    }
    else if( pdgN == -2212){
      fV0tags[iTrackN] = -16;
    }
    
    fV0motherIndex[iTrackN] = iv0;
    fV0motherPDG[iTrackN] = pdgV0;

  }
}


//______________________________________________________________________________
void AliTPCcalibResidualPID::ClearV0PIDlist()
{
  //
  // Clear the PID tag list
  //

  delete [] fV0tags;
  fV0tags = 0;
  
  delete [] fV0motherIndex;
  fV0motherIndex = 0;
  
  delete [] fV0motherPDG;
  fV0motherPDG = 0;
  
  fNumTagsStored = 0;
}


//______________________________________________________________________________
Char_t AliTPCcalibResidualPID::GetV0tag(Int_t trackIndex) const
{
  //
  // Get the tag for the corresponding trackIndex. Returns -99 in case of invalid index/tag list.
  //
  
  if (trackIndex < 0 || trackIndex >= fNumTagsStored || !fV0tags)
    return -99;
  
  return fV0tags[trackIndex];
}


//______________________________________________________________________________
Int_t AliTPCcalibResidualPID::GetV0motherIndex(Int_t trackIndex) const
{
  //
  // Get the index of the V0 mother for the corresponding trackIndex. Returns -99 in case of invalid index/mother index list.
  //
  
  if (trackIndex < 0 || trackIndex >= fNumTagsStored || !fV0motherIndex)
    return -99;
  
  return fV0motherIndex[trackIndex];
}


//______________________________________________________________________________
Int_t AliTPCcalibResidualPID::GetV0motherPDG(Int_t trackIndex) const
{
  //
  // Get the PDG code of the V0 mother for the corresponding trackIndex. Returns 0 in case of invalid index/mother index list.
  //
  
  if (trackIndex < 0 || trackIndex >= fNumTagsStored || !fV0motherPDG)
    return 0;
  
  return fV0motherPDG[trackIndex];
}


//______________________________________________________________________________
Int_t AliTPCcalibResidualPID::MergeGraphErrors(TGraphErrors* mergedGraph, TCollection* li) 
{
  // Adds all graphs with errors from the collection to this graph.
  // Returns the total number of poins in the result or -1 in case of an error.
  
  // TODO "Normal" merge of latest root trunk will also take into account the error bars,
  // so this function can be replaced by the normal root function with the latest root version.
  
  TIter next(li);
  while (TObject* o = next()) {
    TGraph *g = dynamic_cast<TGraph*>(o);
    if (!g) {
      Printf("ERROR: Cannot merge an object which doesn't inherit from TGraph found in the list");
      return -1;
    }
    int n0 = mergedGraph->GetN();
    int n1 = n0+g->GetN();
    mergedGraph->Set(n1);
    Double_t * x = g->GetX();
    Double_t * y = g->GetY();
    Double_t * ex = g->GetEX();
    Double_t * ey = g->GetEY();
    for (Int_t i = 0 ; i < g->GetN(); i++) {
      mergedGraph->SetPoint(n0+i, x[i], y[i]);
      Double_t exPoint = ex ? ex[i] : 0;
      Double_t eyPoint = ey ? ey[i] : 0;
      mergedGraph->SetPointError(n0+i, exPoint, eyPoint);
    }
  }
  return mergedGraph->GetN();
}

//________________________________________________________________________
Bool_t AliTPCcalibResidualPID::TPCCutMIGeo(const AliVTrack* track, const AliVEvent* evt, TTreeStream* streamer)
{
  //
  // TPC Cut MIGeo
  //

  if (!track || !evt)
    return kFALSE;
  
  const Short_t sign = track->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[100];

  track->GetXYZ(xyz);
  track->GetPxPyPz(pxpypz);

  AliExternalTrackParam* par = new AliExternalTrackParam(xyz, pxpypz, cv, sign);
  const AliESDtrack dummy;

  const Double_t magField = evt->GetMagneticField();
  Double_t varGeom = dummy.GetLengthInActiveZone(par, 3, 236, magField, 0, 0);
  Double_t varNcr  = track->GetTPCClusterInfo(3, 1);
  Double_t varNcls = track->GetTPCsignalN();

  const Double_t varEval = 130. - 5. * TMath::Abs(1. / track->Pt());
  Bool_t cutGeom   = varGeom > fgCutGeo * varEval;
  Bool_t cutNcr    = varNcr  > fgCutNcr * varEval;
  Bool_t cutNcls   = varNcls > fgCutNcl * varEval;

  Bool_t kout = cutGeom && cutNcr && cutNcls;

  if (streamer) {
    Double_t dedx = track->GetTPCsignal();
//     Double_t dedx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);

    (*streamer)<<"tree"<<
      "param.="<< par<<
      "varGeom="<<varGeom<<
      "varNcr="<<varNcr<<
      "varNcls="<<varNcls<<
      
      "cutGeom="<<cutGeom<<
      "cutNcr="<<cutNcr<<
      "cutNcls="<<cutNcls<<
      
      "kout="<<kout<<
      "dedx="<<dedx<<
      
      "\n";
  }
  
  delete par;
  
  return kout;
}

//________________________________________________________________________
void AliTPCcalibResidualPID::SetAxisNamesFromTitle(const THnSparseF *h)
{
  // set the histogram axis names from the axis titles
  for (Int_t i=0; i<h->GetNdimensions(); ++i) {
    h->GetAxis(i)->SetName(h->GetAxis(i)->GetTitle());
  }
}

//________________________________________________________________________
TString AliTPCcalibResidualPID::GetStringFitType(Int_t fitType) {   
  if (fitType == AliTPCcalibResidualPID::kLund)
    return "Lund";
  else if (fitType == AliTPCcalibResidualPID::kSaturatedLund)
    return "Saturated Lund";
  else if (fitType == AliTPCcalibResidualPID::kAleph)
    return "ALEPH";
  else if (fitType == AliTPCcalibResidualPID::kAlephWithAdditionalParam)
    return "Modified ALEPH";
  else if (fitType == AliTPCcalibResidualPID::kAlephExternal)
    return "ExternalTrackParameters::ALEPH";
  else
    return "";
}
