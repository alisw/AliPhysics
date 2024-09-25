/**************************************************************************
 * Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

//--------------------------------------------------------------------------------
// PID CME fast analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
// ROOT classes
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TSpline.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskCMEPIDCVE.h"

ClassImp(AliAnalysisTaskCMEPIDCVE);

//---------------------------------------------------
AliAnalysisTaskCMEPIDCVE::AliAnalysisTaskCMEPIDCVE() :
  AliAnalysisTaskSE(),
  fMASS_LAMBDA(1.115683),
  isTightPileUp(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fFilterBit(96),
  fNclsCut(70),
  fChi2Max(2.5),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaCut(0.8),
  fDedxCut(10.),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fAOD(nullptr),
  fPIDResponse(nullptr),
  fUtils(nullptr),
  fMultSel(nullptr),
  fRunNum(-999),
  fOldRunNum(-999),
  fRunNumBin(-999),
  fCent(-999),
  fCentBin(-999),
  fCentV0M(-999),
  fCentTRK(-999),
  fCentSPD0(-999),
  fCentSPD1(-999),
  fSumQ2xTPC(-999),
  fSumQ2yTPC(-999),
  fWgtMultTPC(-999),
  mapTPCTrksIDPhiWgt(0),
  vecParticle(0),
  vecLambda(0),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  listQA(nullptr),
  hEvtCount(nullptr),
  hRunNumBin(nullptr),
  hPt(nullptr),
  hEta(nullptr),
  hPhi(nullptr),
  hNhits(nullptr),
  h2PDedx(nullptr),
  hDcaXY(nullptr),
  hDcaZ(nullptr),
  listResults(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 2; i++) hCent[i] = nullptr;
  for (int i = 0; i < 2; i++) hVz[i] = nullptr;
  for (int i = 0; i < 8; i++) h2Cent[i] = nullptr;
  for (int i = 0; i < 2; i++) h2MultCent[i] = nullptr;
  for (int i = 0; i < 6; i++) h2MultMult[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtPion[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaPion[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiPion[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtProton[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaProton[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiProton[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hRapLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hMassLambda[i] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaOS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaSS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pGammaOS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pGammaSS[i][j] = nullptr;
  pDeltaOS_LambdaProton = nullptr;
  pDeltaSS_LambdaProton = nullptr;
  pDeltaOS_LambdaPion = nullptr;
  pDeltaSS_LambdaPion = nullptr;
  pGammaOS_LambdaProton = nullptr;
  pGammaSS_LambdaProton = nullptr;
  pGammaOS_LambdaPion = nullptr;
  pGammaSS_LambdaPion = nullptr;
  pDeltaOS_LambdaLambda = nullptr;
  pDeltaSS_LambdaLambda = nullptr;
  pGammaOS_LambdaLambda = nullptr;
  pGammaSS_LambdaLambda = nullptr;
  for (int i = 0; i < 3; i++) p3DeltaCentMassMass[i] = nullptr;
  for (int i = 0; i < 3; i++) p3GammaCentMassMass[i] = nullptr;
  for (int i = 0; i < 3; i++) h3LambdaCentMassMass[i] = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskCMEPIDCVE::AliAnalysisTaskCMEPIDCVE(const char *name) :
  AliAnalysisTaskSE(name),
  fMASS_LAMBDA(1.115683),
  isTightPileUp(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fFilterBit(96),
  fNclsCut(70),
  fChi2Max(2.5),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaCut(0.8),
  fDedxCut(10.),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fAOD(nullptr),
  fPIDResponse(nullptr),
  fUtils(nullptr),
  fMultSel(nullptr),
  fRunNum(-999),
  fOldRunNum(-999),
  fRunNumBin(-999),
  fCent(-999),
  fCentBin(-999),
  fCentV0M(-999),
  fCentTRK(-999),
  fCentSPD0(-999),
  fCentSPD1(-999),
  fSumQ2xTPC(-999),
  fSumQ2yTPC(-999),
  fWgtMultTPC(-999),
  mapTPCTrksIDPhiWgt(0),
  vecParticle(0),
  vecLambda(0),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  listQA(nullptr),
  hEvtCount(nullptr),
  hRunNumBin(nullptr),
  hPt(nullptr),
  hEta(nullptr),
  hPhi(nullptr),
  hNhits(nullptr),
  h2PDedx(nullptr),
  hDcaXY(nullptr),
  hDcaZ(nullptr),
  listResults(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 2; i++) hCent[i] = nullptr;
  for (int i = 0; i < 2; i++) hVz[i] = nullptr;
  for (int i = 0; i < 8; i++) h2Cent[i] = nullptr;
  for (int i = 0; i < 2; i++) h2MultCent[i] = nullptr;
  for (int i = 0; i < 6; i++) h2MultMult[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtPion[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaPion[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiPion[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtProton[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaProton[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiProton[i] = nullptr;
  for (int i = 0; i < 2; i++) hPtLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hEtaLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hPhiLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hRapLambda[i] = nullptr;
  for (int i = 0; i < 2; i++) hMassLambda[i] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaOS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaSS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pGammaOS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pGammaSS[i][j] = nullptr;
  pDeltaOS_LambdaProton = nullptr;
  pDeltaSS_LambdaProton = nullptr;
  pDeltaOS_LambdaPion = nullptr;
  pDeltaSS_LambdaPion = nullptr;
  pGammaOS_LambdaProton = nullptr;
  pGammaSS_LambdaProton = nullptr;
  pGammaOS_LambdaPion = nullptr;
  pGammaSS_LambdaPion = nullptr;
  pDeltaOS_LambdaLambda = nullptr;
  pDeltaSS_LambdaLambda = nullptr;
  pGammaOS_LambdaLambda = nullptr;
  pGammaSS_LambdaLambda = nullptr;
  for (int i = 0; i < 3; i++) p3DeltaCentMassMass[i] = nullptr;
  for (int i = 0; i < 3; i++) p3GammaCentMassMass[i] = nullptr;
  for (int i = 0; i < 3; i++) h3LambdaCentMassMass[i] = nullptr;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskCMEPIDCVE::~AliAnalysisTaskCMEPIDCVE()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (listQA) delete listQA;
  if (listResults) delete listResults;
}

//---------------------------------------------------

void AliAnalysisTaskCMEPIDCVE::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskCMEPIDCVE::UserCreateOutputObjects()
{
  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Rihan 18q/r Pile-up function
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

    double parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
    fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    double parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
    fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    double parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
    fMultCutPU->SetParameters(parFB32);
  }

  //------------------
  // QA
  //------------------
  listQA = new TList();
  listQA -> SetName("listQA");
  listQA -> SetOwner(kTRUE);
  // event-wise
  hEvtCount = new TH1D("EvtCount", "Event Count", 20, 1, 21);
  listQA->Add(hEvtCount);

  ////////////////////////
  // Run Number Info
  ////////////////////////
  hRunNumBin = new TH1I("runNum","Counts of Run Number;Run Number;Counts", 1, 0, 1);
  listQA->Add(hRunNumBin);

  // Event-wise QA
  hCent[0] = new TH1D("hCentBfCut", "Dist. of Centrality Before Cut", 80, 0., 80.);
  hCent[1] = new TH1D("hCentAfCut", "Dist. of Centrality After Cut", 80, 0., 80.);

  hVz[0] = new TH1D("hVzBfCut", "Dist of Vz Before Cut", 200, -10., 10.);
  hVz[1] = new TH1D("hVzAfCut", "Dist of Vz After Cut", 200, -10., 10.);

  h2Cent[0] = new TH2D("h2Cent_V0M_SPD1_BfCut", ";centV0M;centSPD1", 80, 0, 80., 80, 0, 80.);
  h2Cent[1] = new TH2D("h2Cent_V0M_SPD1_AfCut", ";centV0M;centSPD1", 80, 0, 80., 80, 0, 80.);
  h2Cent[2] = new TH2D("h2Cent_V0M_TRK_BfCut", ";centV0M;centTRK", 80, 0, 80., 80, 0, 80.);
  h2Cent[3] = new TH2D("h2Cent_V0M_TRK_AfCut", ";centV0M;centTRK", 80, 0, 80., 80, 0, 80.);
  h2Cent[4] = new TH2D("h2Cent_V0M_SPD0_BfCut", ";centV0M;centSPD0", 80, 0, 80., 80, 0, 80.);
  h2Cent[5] = new TH2D("h2Cent_V0M_SPD0_AfCut", ";centV0M;centSPD0", 80, 0, 80., 80, 0, 80.);
  h2Cent[6] = new TH2D("h2Cent_SPD1_SPD0_BfCut", ";centSPD1;centSPD0", 80, 0, 80., 80, 0, 80.);
  h2Cent[7] = new TH2D("h2Cent_SPD1_SPD0_AfCut", ";centSPD1;centSPD0", 80, 0, 80., 80, 0, 80.);

  h2MultCent[0] = new TH2D("h2MultCent_BfCut", ";centV0M;multFB32", 80, 0, 80., 50, 0, 5000.);
  h2MultCent[1] = new TH2D("h2MultCent_AfCut", ";centV0M;multFB32", 80, 0, 80., 50, 0, 5000.);

  listQA->Add(hCent[0]);
  listQA->Add(hCent[1]);
  listQA->Add(hVz[0]);
  listQA->Add(hVz[1]);
  for (int i = 0; i < 8; i++) listQA->Add(h2Cent[i]);
  listQA->Add(h2MultCent[0]);
  listQA->Add(h2MultCent[1]);

  // track-wise QA
  hPt  = new TH1D("hPt", ";p_{T}", 200, 0., 20.);
  hEta = new TH1D("hEta", ";#eta", 100, -2., 2.);
  hPhi = new TH1D("hPhi", ";#phi", 100, 0, TMath::TwoPi());
  hNhits = new TH1D("hNhits", ";nhits", 200, 0., 200.);
  h2PDedx = new TH2D("h2PDedx", ";pdedx", 400, -10., 10., 400, 0, 1000);
  hDcaXY = new TH1D("hDcaXY", ";DcaXY", 500, 0., 1);
  hDcaZ  = new TH1D("hDcaZ",  ";DcaZ", 500, 0., 1);
  listQA->Add(hPt);
  listQA->Add(hEta);
  listQA->Add(hPhi);
  listQA->Add(hNhits);
  listQA->Add(h2PDedx);
  listQA->Add(hDcaXY);
  listQA->Add(hDcaZ);


  for (int i = 0; i < 2; i++) hPtPion[i] = new TH1D(Form("hPtPion_%i",i), ";p_{T}", 200, 0., 20.);
  for (int i = 0; i < 2; i++) hEtaPion[i] = new TH1D(Form("hEtaPion_%i",i), ";#eta", 100, -2., 2.);
  for (int i = 0; i < 2; i++) hPhiPion[i] = new TH1D(Form("hPhiPion_%i",i), ";#phi", 100, 0, TMath::TwoPi());

  for (int i = 0; i < 2; i++) hPtKaon[i] = new TH1D(Form("hPtKaon_%i",i), ";p_{T}", 200, 0., 20.);
  for (int i = 0; i < 2; i++) hEtaKaon[i] = new TH1D(Form("hEtaKaon_%i",i), ";#eta", 100, -2., 2.);
  for (int i = 0; i < 2; i++) hPhiKaon[i] = new TH1D(Form("hPhiKaon_%i",i), ";#phi", 100, 0, TMath::TwoPi());

  for (int i = 0; i < 2; i++) hPtProton[i] = new TH1D(Form("hPtProton_%i",i), ";p_{T}", 200, 0., 20.);
  for (int i = 0; i < 2; i++) hEtaProton[i] = new TH1D(Form("hEtaProton_%i",i), ";#eta", 100, -2., 2.);
  for (int i = 0; i < 2; i++) hPhiProton[i] = new TH1D(Form("hPhiProton_%i",i), ";#phi", 100, 0, TMath::TwoPi());

  for (int i = 0; i < 2; i++) hPtLambda[i] = new TH1D(Form("hPtLambda_%i",i), ";p_{T}", 200, 0., 20.);
  for (int i = 0; i < 2; i++) hEtaLambda[i] = new TH1D(Form("hEtaLambda_%i",i), ";#eta", 100, -2., 2.);
  for (int i = 0; i < 2; i++) hPhiLambda[i] = new TH1D(Form("hPhiLambda_%i",i), ";#phi", 100, 0, TMath::TwoPi());
  for (int i = 0; i < 2; i++) hRapLambda[i] = new TH1D(Form("hRapLambda_%i",i), ";y", 100, -2., 2.);
  for (int i = 0; i < 2; i++) hMassLambda[i] = new TH1D(Form("hMassLambda_%i",i), ";mass", 100, 1.11, 1.12);

  for (int i = 0; i < 2; i++) listQA->Add(hPtPion[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaPion[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiPion[i]);

  for (int i = 0; i < 2; i++) listQA->Add(hPtKaon[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaKaon[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiKaon[i]);

  for (int i = 0; i < 2; i++) listQA->Add(hPtProton[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaProton[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiProton[i]);

  for (int i = 0; i < 2; i++) listQA->Add(hPtLambda[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaLambda[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiLambda[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hRapLambda[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hMassLambda[i]);

  PostData(1,listQA);
  if (fDebug) Printf("Post fQAList Data Success!");

  ////////////////////////
  // Results
  ////////////////////////
  listResults = new TList();
  listResults -> SetName("listResults");
  listResults -> SetOwner(kTRUE);
  //        pion kaon proton
  //pion    00   01   02
  //kaon    10   11   12
  //proton  20   21   22

  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaOS[i][j] = new TProfile(Form("pDeltaOS_%i%i",i,j), ";centrality;#delta_{OS}", 7, 0, 70);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaSS[i][j] = new TProfile(Form("pDeltaSS_%i%i",i,j), ";centrality;#delta_{SS}", 7, 0, 70);
  // for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) listResults->Add(pDeltaOS[i][j]);
  // for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) listResults->Add(pDeltaSS[i][j]);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pGammaOS[i][j] = new TProfile(Form("pGammaOS_%i%i",i,j), ";centrality;#gamma_{OS}", 7, 0, 70);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pGammaSS[i][j] = new TProfile(Form("pGammaSS_%i%i",i,j), ";centrality;#gamma_{SS}", 7, 0, 70);
  // for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) listResults->Add(pGammaOS[i][j]);
  // for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) listResults->Add(pGammaSS[i][j]);

  pDeltaOS_LambdaProton = new TProfile("pDeltaOS_LambdaProton", ";centrality;#delta_{OS}", 7, 0, 70);
  pDeltaSS_LambdaProton = new TProfile("pDeltaSS_LambdaProton", ";centrality;#delta_{SS}", 7, 0, 70);
  pDeltaOS_LambdaPion = new TProfile("pDeltaOS_LambdaPion", ";centrality;#delta_{OS}", 7, 0, 70);
  pDeltaSS_LambdaPion = new TProfile("pDeltaSS_LambdaPion", ";centrality;#delta_{SS}", 7, 0, 70);
  listResults->Add(pDeltaOS_LambdaProton);
  listResults->Add(pDeltaSS_LambdaProton);
  listResults->Add(pDeltaOS_LambdaPion);
  listResults->Add(pDeltaSS_LambdaPion);

  pGammaOS_LambdaProton = new TProfile("pGammaOS_LambdaProton", ";centrality;#gamma_{OS}", 7, 0, 70);
  pGammaSS_LambdaProton = new TProfile("pGammaSS_LambdaProton", ";centrality;#gamma_{SS}", 7, 0, 70);
  pGammaOS_LambdaPion = new TProfile("pGammaOS_LambdaPion", ";centrality;#gamma_{OS}", 7, 0, 70);
  pGammaSS_LambdaPion = new TProfile("pGammaSS_LambdaPion", ";centrality;#gamma_{SS}", 7, 0, 70);
  listResults->Add(pGammaOS_LambdaProton);
  listResults->Add(pGammaSS_LambdaProton);
  listResults->Add(pGammaOS_LambdaPion);
  listResults->Add(pGammaSS_LambdaPion);

  pDeltaOS_LambdaLambda = new TProfile("pDeltaOS_LambdaLambda", ";centrality;#delta_{OS}", 7, 0, 70);
  pDeltaSS_LambdaLambda = new TProfile("pDeltaSS_LambdaLambda", ";centrality;#delta_{SS}", 7, 0, 70);
  pGammaOS_LambdaLambda = new TProfile("pGammaOS_LambdaLambda", ";centrality;#gamma_{OS}", 7, 0, 70);
  pGammaSS_LambdaLambda = new TProfile("pGammaSS_LambdaLambda", ";centrality;#gamma_{SS}", 7, 0, 70);
  listResults->Add(pDeltaOS_LambdaLambda);
  listResults->Add(pDeltaSS_LambdaLambda);
  listResults->Add(pGammaOS_LambdaLambda);
  listResults->Add(pGammaSS_LambdaLambda);

  //obv - (cent - mass - mass)
  for (int i = 0; i < 3; i++) p3DeltaCentMassMass[i] = new TProfile3D(Form("p3DeltaCentMassMass_%i",i), ";centrality;mass_{p#pi};mass_{p#pi}", 7, 0, 70, 15, fMASS_LAMBDA - 0.02, fMASS_LAMBDA + 0.02, 15, fMASS_LAMBDA - 0.02, fMASS_LAMBDA + 0.02);
  for (int i = 0; i < 3; i++) p3GammaCentMassMass[i] = new TProfile3D(Form("p3GammaCentMassMass_%i",i), ";centrality;mass_{p#pi};mass_{p#pi}", 7, 0, 70, 15, fMASS_LAMBDA - 0.02, fMASS_LAMBDA + 0.02, 15, fMASS_LAMBDA - 0.02, fMASS_LAMBDA + 0.02);
  for (int i = 0; i < 3; i++) h3LambdaCentMassMass[i] = new TH3D(Form("h3LambdaCentMassMass_%i",i), ";centrality;mass_{p#pi};mass_{#pi#pi}", 7, 0, 70, 15, fMASS_LAMBDA - 0.02, fMASS_LAMBDA + 0.02, 15, fMASS_LAMBDA - 0.02, fMASS_LAMBDA + 0.02);
  for (int i = 0; i < 3; i++) listResults->Add(p3DeltaCentMassMass[i]);
  for (int i = 0; i < 3; i++) listResults->Add(p3GammaCentMassMass[i]);
  for (int i = 0; i < 3; i++) listResults->Add(h3LambdaCentMassMass[i]);

  PostData(2,listResults);
  if (fDebug) Printf("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskCMEPIDCVE::UserExec(Option_t *)
{
  if (fDebug) Printf("===============================We are in UserExec!!!===================================");
  hEvtCount->Fill("All",1);
  //----------------------------
  // Handle
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else hEvtCount->Fill("Manager",1);

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else hEvtCount->Fill("Handler",1);

  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else hEvtCount->Fill("fAOD",1);

  fPIDResponse = handler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else hEvtCount->Fill("fPID",1);

  fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else hEvtCount->Fill("fUtils",1);

  fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  if (!fMultSel) {
    AliError(Form("%s: Could not get AliMultSelection", GetName()));
  } else hEvtCount->Fill("fMultSel",1);

  if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils || !fMultSel) return;
  if (fDebug) Printf("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  unsigned int mask = handler->IsEventSelected();
  //bool isTrigselected = false;
  // if (fTrigger.EqualTo("kMB"))
  // isTrigselected = mask & AliVEvent::kMB;
  // else if (fTrigger.EqualTo("kINT7"))
  // isTrigselected = mask & AliVEvent::kINT7;
  // else if (fTrigger.EqualTo("kINT7+kSemiCentral"))
  // isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral);
  // else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
  // isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  // if (isTrigselected == false) return;
  if (!(mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral))) return;
  hEvtCount->Fill("Trigger",1);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Run Number
  //----------------------------
  fRunNum = fAOD->GetRunNumber();
  hEvtCount->Fill("Run Number",1);
  hRunNumBin->Fill(TString(fRunNum),1);
  if (fDebug) Printf("read in done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  fVtx -> GetXYZ(fVertex);
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (fabs(fVertex[0])<1e-6 || fabs(fVertex[1])<1e-6 || fabs(fVertex[2])<1e-6) return;
  hVz[0]->Fill(fVertex[2]);
  if (fabs(fVertex[2]) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  hVz[1]->Fill(fVertex[2]);
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  fCentV0M  = fMultSel->GetMultiplicityPercentile("V0M");
  fCentTRK  = fMultSel->GetMultiplicityPercentile("TRK");
  fCentSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
  fCentSPD1 = fMultSel->GetMultiplicityPercentile("CL1");

  //we use centV0M as the default centrality
  fCent = fCentV0M;
  h2Cent[0]->Fill(fCentV0M,fCentSPD1);
  h2Cent[2]->Fill(fCentV0M,fCentTRK);
  h2Cent[4]->Fill(fCentV0M,fCentSPD0);
  h2Cent[6]->Fill(fCentSPD1,fCentSPD0);
  if (fabs(fCent-fCentSPD1)>fCentDiffCut) return;
  h2Cent[1]->Fill(fCentV0M,fCentSPD1);
  h2Cent[3]->Fill(fCentV0M,fCentTRK);
  h2Cent[5]->Fill(fCentV0M,fCentSPD0);
  h2Cent[7]->Fill(fCentSPD1,fCentSPD0);
  if (fCent < 0 || fCent >= 80) return;

  // cent
  hCent[0]->Fill(fCent);
  hEvtCount->Fill("Centrality",1);
  if (fDebug) Printf("centrality done!");

  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    if (!RejectEvtTFFit()) return;
  }
  hCent[1]->Fill(fCent);
  hEvtCount->Fill("Pile-up",1);
  if (fDebug) Printf("pile-up done!");
  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  //Reset vectors
  ResetVectors();
  hEvtCount->Fill("Reset Vectors",1);
  // must loop tracks becasue we need the TPC plane
  if (!LoopTracks()) return;
  hEvtCount->Fill("Loop Tracks",1);
  if (fDebug) Printf("Loop Tracks done!");
  //----------------------------
  // Loop V0s
  //----------------------------
  if (!LoopV0s()) return;
  hEvtCount->Fill("Loop V0s",1);
  if (fDebug) Printf("Loop V0s done!");
  //----------------------------
  // Pair
  //----------------------------
  // if (!PairTrkTrk()) return;
  // hEvtCount->Fill("Pair Trk-Trk",1);
  // if (fDebug) Printf("Pair V0 & Trk done!");
  // if (!PairV0Trk()) return;
  // hEvtCount->Fill("Pair V0-Trk",1);
  // if (fDebug) Printf("Pair V0 & Trk done!");
  if (!PairV0V0()) return;
  hEvtCount->Fill("Pair V0-V0",1);
  if (fDebug) Printf("Pair V0 & V0 done!");
  //------------------
  // Post output data.
  //------------------
  PostData(1,listQA);
  PostData(2,listResults);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::LoopTracks()
{
  int nTrks = fAOD->GetNumberOfTracks();
  if (nTrks < 10) return false;
  for (int iTrk = 0; iTrk < nTrks; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (!AcceptAODTrack(track)) continue;
    //------------------
    // NUE & NUA
    //------------------
    double  phi = track->Phi();
    double   pt = track->Pt();
    double  eta = track->Eta();
    int  charge = track->Charge();
    int   nhits = track->GetTPCNcls();
    double dedx = track->GetTPCsignal();
    int      id = track->GetID();
    
    hPt->Fill(pt);
    hEta->Fill(eta);
    hPhi->Fill(phi);
    hNhits->Fill(nhits);
    h2PDedx->Fill(track->P()*charge, dedx);

    //DCA Cut
    double dcaxy = -999, dcaz = -999;
    if(!GetDCA(dcaxy,dcaz,track)) continue;
    //just use FB96
    // if (fabs(dcaxy) > fDcaCutXY) continue;
    // if (fabs(dcaz) > fDcaCutZ) continue;
    hDcaXY->Fill(fabs(dcaxy));
    hDcaZ->Fill(fabs(dcaz));

    // TPC Plane
    double weight = 1.0;
    if (pt > 0.2 && pt < 2.0) {
      fSumQ2xTPC  += weight * cos(2 * phi);
      fSumQ2yTPC  += weight * sin(2 * phi);
      fWgtMultTPC += weight;
      std::vector<double> vec_phi_weight;
      vec_phi_weight.emplace_back(phi);
      vec_phi_weight.emplace_back(weight);
      mapTPCTrksIDPhiWgt[id] = vec_phi_weight;
    }

    int pid = GetPIDofParticle(track);
    if (std::isnan(pid)) continue;
    int bMatter = 0;
    if (pid < 0) bMatter = 1;

    if (abs(pid) == 211) {
      hPtPion[bMatter] -> Fill(pt);
      hEtaPion[bMatter] -> Fill(eta);
      hPhiPion[bMatter] -> Fill(phi);
    }
    else if (abs(pid) == 321) {
      hPtKaon[bMatter] -> Fill(pt);
      hEtaKaon[bMatter] -> Fill(eta);
      hPhiKaon[bMatter] -> Fill(phi);
    }
    else if (abs(pid) == 2212) {
      hPtProton[bMatter] -> Fill(pt);
      hEtaProton[bMatter] -> Fill(eta);
      hPhiProton[bMatter] -> Fill(phi);
    }
    else continue;

    vecParticle.emplace_back(track,pid);
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::PairTrkTrk()
{
  for (auto particle_1 : vecParticle) {
    double pt_1 = particle_1.first->Pt();
    double eta_1 = particle_1.first->Eta();
    double phi_1 = particle_1.first->Phi();
    int id_1 = particle_1.first->GetID();
    int pid_1 = particle_1.second;

    for (auto particle_2 : vecParticle) {
      double pt_2 = particle_2.first->Pt();
      double eta_2 = particle_2.first->Eta();
      double phi_2 = particle_2.first->Phi();
      int id_2 = particle_2.first->GetID();
      int pid_2 = particle_2.second;

      if(id_1 == id_2) continue;
      double delta = cos(phi_1 - phi_2);
      double psi = GetTPCPlaneNoAutoCorr({id_1, id_2});
      if (std::isnan(psi)) continue;
      double gamma = cos(phi_1 + phi_2 - 2 * psi);

      //pid   211  321 2212
      //211    00   01   02
      //321    10   11   12
      //2212   20   21   22

      int row_pid = -999;
      int col_pid = -999;
      if(abs(pid_1) == 211 && abs(pid_2) == 211)        {row_pid = 0; col_pid = 0;}
      else if(abs(pid_1) ==  211 && abs(pid_2) ==  321) {row_pid = 0; col_pid = 1;}
      else if(abs(pid_1) ==  211 && abs(pid_2) == 2212) {row_pid = 0; col_pid = 2;}
      else if(abs(pid_1) ==  321 && abs(pid_2) ==  211) {row_pid = 1; col_pid = 0;}
      else if(abs(pid_1) ==  321 && abs(pid_2) ==  321) {row_pid = 1; col_pid = 1;}
      else if(abs(pid_1) ==  321 && abs(pid_2) == 2212) {row_pid = 1; col_pid = 2;}
      else if(abs(pid_1) == 2212 && abs(pid_2) ==  211) {row_pid = 2; col_pid = 0;}
      else if(abs(pid_1) == 2212 && abs(pid_2) ==  321) {row_pid = 2; col_pid = 1;}
      else if(abs(pid_1) == 2212 && abs(pid_2) == 2212) {row_pid = 2; col_pid = 2;}
      else continue;

      if(pid_1 * pid_2 < 0) {
        pDeltaOS[row_pid][col_pid]->Fill(fCent,delta);
        pGammaOS[row_pid][col_pid]->Fill(fCent,gamma);
      } else if(pid_1 * pid_2 > 0) {
        pDeltaSS[row_pid][col_pid]->Fill(fCent,delta);
        pGammaSS[row_pid][col_pid]->Fill(fCent,gamma);
      } else continue;
    }
  }
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskCMEPIDCVE::ResetVectors()
{
  fSumQ2xTPC = 0;
  fSumQ2yTPC = 0;
  fWgtMultTPC = 0;
  mapTPCTrksIDPhiWgt.clear();
  std::unordered_map<int, std::vector<double>>().swap(mapTPCTrksIDPhiWgt);
  vecParticle.clear();
  std::vector<std::pair<AliAODTrack*,int>>().swap(vecParticle);
  vecLambda.clear();
  std::vector<std::pair<AliAODv0*,int>>().swap(vecLambda);
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::RejectEvtTFFit()
{
  int nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  int nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  int nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  int nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const int nTracks = fAOD->GetNumberOfTracks();
  int multTrk = 0;
  for (int it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk) {
          delete aodTrk;
          continue;
      }
      if (aodTrk->TestFilterBit(32)) {
        if ((abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  h2MultCent[0]->Fill(fCentV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  float  multV0a = aodV0->GetMTotV0A();
  float  multV0c = aodV0->GetMTotV0C();
  float  multV0Tot = multV0a + multV0c;
  unsigned short multV0aOn = aodV0->GetTriggerChargeA();
  unsigned short multV0cOn = aodV0->GetTriggerChargeC();
  unsigned short multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (fCentSPD0 < fCenCutLowPU->Eval(fCentV0M)) return false;
  if (fCentSPD0 > fCenCutHighPU->Eval(fCentV0M)) return false;
  if (float(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (float(multTrk) < fMultCutPU->Eval(fCentV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  if (isTightPileUp == true){
    int tpcClsTot = fAOD->GetNumberOfTPCClusters();
    float nclsDif = float(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot);
    if (nclsDif > 200000)//can be varied to 150000, 200000
    return false;
  }

  h2MultCent[1]->Fill(fCentV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}


//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::AcceptAODTrack(AliAODTrack *track)
{
    //------------------
    // track cut
    //------------------
    double pt = track->Pt();
    if(pt < fPtMin) return false;
    if(pt > fPtMax) return false;

    double eta = track->Eta();
    if(fabs(eta)  > fEtaCut) return false;

    int nhits  = track->GetTPCNcls();
    if(fabs(nhits) < fNclsCut) return false;

    double dedx = track->GetTPCsignal();
    if(dedx < fDedxCut) return false;

    double chi2   = track->Chi2perNDF();
    if(chi2 < fChi2Min) return false;
    if(chi2 > fChi2Max) return false;

    return true;
}

//---------------------------------------------------
int AliAnalysisTaskCMEPIDCVE::GetPIDofParticle(AliAODTrack* ftrack)
{
  if(!fPIDResponse) {
    std::cout<<"AliPIDResponse does not exist!!!"<<std::endl;
    return nan("");
  }

  short charge = ftrack->Charge();
  if (abs(charge) != 1) return 0;
  float nSigTPC = -999;
  float nSigTOF = -999;
  float nSigRMS = -999;

  bool isPion = false;
  bool isKaon = false;
  bool isProton = false;

  //first proton
  nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
  nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
  nSigRMS = hypot(nSigTPC, nSigTOF);
  if (TMath::Abs(nSigRMS) < fNSigmaTOFCut) isProton = true;

  //second kaon
  nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
  nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
  nSigRMS = hypot(nSigTPC, nSigTOF);

  if (ftrack->Pt() <= 0.45 && TMath::Abs(nSigTPC) < fNSigmaTPCCut) isKaon = true;
  if (ftrack->Pt() >  0.45 && TMath::Abs(nSigRMS) < fNSigmaTOFCut) isKaon = true;

  //third pion
  nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
  nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
  nSigRMS = hypot(nSigTPC, nSigTOF);

  if (ftrack->Pt() <= 0.5 && TMath::Abs(nSigTPC) < fNSigmaTPCCut) isPion = true;
  if (ftrack->Pt() >  0.5 && TMath::Abs(nSigRMS) < fNSigmaTOFCut) isPion = true;

  if ((isProton + isKaon + isPion) != 1) return nan("");
  if (isPion)   return 211 * charge;
  if (isProton) return 2212 * charge;
  if (isKaon)   return 321 * charge;
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track)
{
  if(!track) return false;
  double r[3];
  if (track->GetXYZ(r)) {
    dcaxy = r[0];
    dcaz  = r[1];
  } else {
    double dcax = r[0] - fVertex[0];
    double dcay = r[1] - fVertex[1];
    dcaz  = r[2] - fVertex[2];
    dcaxy = sqrt(dcax*dcax + dcay*dcay);
  }
  return true;
}

//---------------------------------------------------
  
bool AliAnalysisTaskCMEPIDCVE::LoopV0s()
{
  int nV0s = fAOD->GetNumberOfV0s();
  if (nV0s < 1) return false;
  for (int iV0 = 0; iV0 < nV0s; iV0++) {
    AliAODv0 *v0 = fAOD->GetV0(iV0);
    if (!v0) {
      AliError(Form("%s: Could not get v0s", GetName()));
      continue;
    }
    //Basic kinematic variable
    double pt = v0->Pt();
    if (pt < 0.5) continue;
    double eta = v0->PseudoRapV0();
    if (fabs(eta) > 0.8) continue;

    //V0 cut
    if (!IsGoodV0(v0)) continue;
    //V0 daughters cut
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    if (!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;
    float pDcaPV = v0->DcaPosToPrimVertex();
    float nDcaPV = v0->DcaNegToPrimVertex();
    if ( nDcaPV<0.05 || pDcaPV<0.05) continue;

    int pid = GetLambdaPID(pTrack,nTrack);
    if (abs(pid) != 3122) continue;

    double rap = v0->RapLambda();
    if(rap < -0.5 || rap > 0.5) continue;

    double mass = -999;
    if (pid > 0) mass = v0->MassLambda();
    else mass = v0->MassAntiLambda();
    if( mass < fMASS_LAMBDA - 0.02 || mass > fMASS_LAMBDA + 0.02) continue;
  
    double phi = v0->Phi();
    //QA
    if(pid == 3122) {
      hPtLambda[0]->Fill(pt);
      hEtaLambda[0]->Fill(eta);
      hPhiLambda[0]->Fill(phi);
      hMassLambda[0]->Fill(mass);
      hRapLambda[0]->Fill(rap);
    } else {
      hPtLambda[1]->Fill(pt);
      hEtaLambda[1]->Fill(eta);
      hPhiLambda[1]->Fill(phi);
      hMassLambda[1]->Fill(mass);
      hRapLambda[1]->Fill(rap);
    }
    vecLambda.emplace_back(v0,pid);
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::IsGoodV0(const AliAODv0 *aodV0)
{
  // Offline reconstructed V0 only
  if (aodV0->GetOnFlyStatus()) return false;
  // Get daughters and check them
  AliAODTrack* myTrackPosTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
  AliAODTrack* myTrackNegTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return false;
  }
  // Unlike signs of daughters
  if (myTrackNegTest->Charge() * myTrackPosTest->Charge() > 0) return false;
  // Cosinus of pointing angle < 0.997
  double dCPA = aodV0->CosPointingAngle(fVertex);
  if (dCPA < 0.997) return false;
  // DCA of V0 < 1.5 cm
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if (abs(dV0Dca) > 1.5) return false;
  // V0 path length before decay 3-100 cm
  double dDecayLength = aodV0->DecayLengthV0(fVertex);
  if (dDecayLength > 100.) return false;
  if (dDecayLength < 3.) return false;
  // DCA between daughters < 0.5cm
  double dDCA = aodV0->DcaV0Daughters();
  if (dDCA > 0.5) return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::IsGoodDaughterTrack(const AliAODTrack *track)
{
  // TPC refit
  if (!track->IsOn(AliAODTrack::kTPCrefit)) return false;
  // No kinks
  if (int(track->GetProdVertex()->GetType()) == AliAODVertex::kKink) return false;
  // Maximum value of transverse momentum
  double dPt = track->Pt();
  if (dPt > 20.) return false;
  // Maximum value of pseudorapidity
  double dEta = track->Eta();
  if (abs(dEta) > 0.8) return false;
  // Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return false;
  // Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return false;
  // [number of crossed rows]>0.8  [number of findable clusters].
  if (nCrossedRowsTPC/findable < 0.8) return false;
  return true;
}

//---------------------------------------------------

int AliAnalysisTaskCMEPIDCVE::GetLambdaPID(const AliAODTrack *pTrack, const AliAODTrack *nTrack)
{
  bool isLambda     = kFALSE;
  bool isAntiLambda = kFALSE;
  int  code = 0;

  //Λ-->(p+)+(π-)
  float nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton));//TPC p+
  float nSigTPCNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));//TPC π-
  //(Λ-)-->(p-)+(π+)
  float nSigTPCPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));//TPC π+
  float nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton));//TPC p-

  isLambda     = (nSigTPCPosProton < 3.) && (nSigTPCNegPion < 3.);
  isAntiLambda = (nSigTPCNegProton < 3.) && (nSigTPCPosPion < 3.);

  if (isLambda)     code =  3122;
  if (isAntiLambda) code = -3122;
  if (isLambda && isAntiLambda) code = 0;
  return code;
}

//---------------------------------------------------
bool AliAnalysisTaskCMEPIDCVE::PairV0Trk() {
  for (auto lambda : vecLambda) {
    double pt_1 = lambda.first->Pt();
    double eta_1 = lambda.first->Eta();
    double phi_1 = lambda.first->Phi();
    int pid_1 = lambda.second;

    int id_daughter_1 = lambda.first->GetPosID();
    int id_daughter_2 = lambda.first->GetNegID();

    for (auto particle : vecParticle) {
      double pt_2 = particle.first->Pt();
      double eta_2 = particle.first->Eta();
      double phi_2 = particle.first->Phi();
      int pid_2 = particle.second;

      int id_2 = particle.first->GetID();
      if (id_2 == id_daughter_1 || id_2 == id_daughter_2) continue;

      double delta = cos(phi_1 - phi_2);
      double psi = GetTPCPlaneNoAutoCorr({id_daughter_1, id_daughter_2, id_2});
      if (std::isnan(psi)) continue;
      double gamma = cos(phi_1 + phi_2 - 2 * psi);

      if(abs(pid_2) == 2212) {
        if (pid_2 * pid_1 < 0)      {
          pDeltaOS_LambdaProton->Fill(fCent,delta);
          pGammaOS_LambdaProton->Fill(fCent,gamma);
        } else if (pid_2 * pid_1 > 0) {
          pDeltaSS_LambdaProton->Fill(fCent,delta);
          pGammaSS_LambdaProton->Fill(fCent,gamma);
        } else continue;
      }
      if(abs(pid_2) == 211) {
        if (pid_2 * pid_1 < 0)      {
          pDeltaOS_LambdaPion->Fill(fCent,delta);
          pGammaOS_LambdaPion->Fill(fCent,gamma);
        } else if (pid_2 * pid_1 > 0) {
          pDeltaSS_LambdaPion->Fill(fCent,delta);
          pGammaSS_LambdaPion->Fill(fCent,gamma);
        } else continue;
      }
    }
  }
  return true;
}

//---------------------------------------------------
double AliAnalysisTaskCMEPIDCVE::GetTPCPlaneNoAutoCorr(std::vector<int> vec_id) {
  double psiNoAuto = nan("");

  double tempSumQ2x = fSumQ2xTPC;
  double tempSumQ2y = fSumQ2yTPC;
  double tempWgtMult = fWgtMultTPC;
  double repeQ2x = 0., repeQ2y = 0., repeWgtMult = 0.;

  std::vector<std::unordered_map<int, std::vector<double>> ::iterator> vec_it;
  for (int id : vec_id) {
    vec_it.push_back(mapTPCTrksIDPhiWgt.find(id));
  }

  for (auto it : vec_it) {
    if(it != mapTPCTrksIDPhiWgt.end()) {
      repeQ2x     += it->second[1] * cos(2 * (it->second)[0]);
      repeQ2y     += it->second[1] * sin(2 * (it->second)[0]);
      repeWgtMult += it->second[1];
    }
  }

  tempSumQ2x  -= repeQ2x;
  tempSumQ2y  -= repeQ2y;
  tempWgtMult -= repeWgtMult;

  if (tempWgtMult > 1.e-6) {
    tempSumQ2x /= tempWgtMult;
    tempSumQ2y /= tempWgtMult;
    psiNoAuto = GetEventPlane(tempSumQ2x,tempSumQ2y,2.);
  } else psiNoAuto = nan("");

  return psiNoAuto;
}

//---------------------------------------------------

inline double AliAnalysisTaskCMEPIDCVE::GetEventPlane(double qx, double qy, double harmonic)
{
  double psi = (1./harmonic)*TMath::ATan2(qy,qx);
  if (psi < 0) return psi += TMath::TwoPi()/harmonic;
  else return psi;
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::PairV0V0()
{
  for (auto lambda_1 : vecLambda) {
    double pt_1 = lambda_1.first->Pt();
    double eta_1 = lambda_1.first->Eta();
    double phi_1 = lambda_1.first->Phi();
    int pid_1 = lambda_1.second;

    int id_daughter_1 = lambda_1.first->GetPosID();
    int id_daughter_2 = lambda_1.first->GetNegID();

    for (auto lambda_2 : vecLambda) {
      double pt_2 = lambda_2.first->Pt();
      double eta_2 = lambda_2.first->Eta();
      double phi_2 = lambda_2.first->Phi();
      int pid_2 = lambda_2.second;

      int id_daughter_3 = lambda_2.first->GetPosID();
      int id_daughter_4 = lambda_2.first->GetNegID();

      if (id_daughter_1 == id_daughter_3 || id_daughter_1 == id_daughter_4) continue;
      if (id_daughter_2 == id_daughter_3 || id_daughter_2 == id_daughter_4) continue;

      double delta = cos(phi_1 - phi_2);
      double psi = GetTPCPlaneNoAutoCorr({id_daughter_1, id_daughter_2, id_daughter_3, id_daughter_4});
      if (std::isnan(psi)) continue;
      double gamma = cos(phi_1 + phi_2 - 2 * psi);

      if (pid_2 * pid_1 < 0) {
        pDeltaOS_LambdaLambda->Fill(fCent,delta);
        pGammaOS_LambdaLambda->Fill(fCent,gamma);
      } else if (pid_2 * pid_1 > 0) {
        pDeltaSS_LambdaLambda->Fill(fCent,delta);
        pGammaSS_LambdaLambda->Fill(fCent,gamma);
      } else continue;

      if (pid_1 == 3122 && pid_2 == 3122) {
        p3DeltaCentMassMass[0] -> Fill(fCent,lambda_1.first->MassLambda(),lambda_2.first->MassLambda(),delta);
        p3GammaCentMassMass[0] -> Fill(fCent,lambda_1.first->MassLambda(),lambda_2.first->MassLambda(),gamma);
        h3LambdaCentMassMass[0] -> Fill(fCent,lambda_1.first->MassLambda(),lambda_2.first->MassLambda());
      } else if (pid_1 == 3122 && pid_2 == -3122) {
        p3DeltaCentMassMass[1] -> Fill(fCent,lambda_1.first->MassLambda(),lambda_2.first->MassAntiLambda(),delta);
        p3GammaCentMassMass[1] -> Fill(fCent,lambda_1.first->MassLambda(),lambda_2.first->MassAntiLambda(),gamma);
        h3LambdaCentMassMass[1] -> Fill(fCent,lambda_1.first->MassLambda(),lambda_2.first->MassAntiLambda());
      } else if (pid_1 == -3122 && pid_2 == -3122) {
        p3DeltaCentMassMass[2] -> Fill(fCent,lambda_1.first->MassAntiLambda(),lambda_2.first->MassAntiLambda(),delta);
        p3GammaCentMassMass[2] -> Fill(fCent,lambda_1.first->MassAntiLambda(),lambda_2.first->MassAntiLambda(),gamma);
        h3LambdaCentMassMass[2] -> Fill(fCent,lambda_1.first->MassAntiLambda(),lambda_2.first->MassAntiLambda());
      }
    }
  }
  return true;
}
