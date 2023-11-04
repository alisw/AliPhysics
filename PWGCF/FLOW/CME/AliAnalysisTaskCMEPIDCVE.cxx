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
  isTightPileUp(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fFilterBit(768),
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
  vecParticle(0),
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
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaOS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaSS[i][j] = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskCMEPIDCVE::AliAnalysisTaskCMEPIDCVE(const char *name) :
  AliAnalysisTaskSE(name),
  isTightPileUp(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fFilterBit(768),
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
  vecParticle(0),
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
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaOS[i][j] = nullptr;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) pDeltaSS[i][j] = nullptr;

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

    Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
    fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
    fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
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

  for (int i = 0; i < 2; i++) listQA->Add(hPtPion[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaPion[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiPion[i]);

  for (int i = 0; i < 2; i++) listQA->Add(hPtKaon[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaKaon[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiKaon[i]);

  for (int i = 0; i < 2; i++) listQA->Add(hPtProton[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hEtaProton[i]);
  for (int i = 0; i < 2; i++) listQA->Add(hPhiProton[i]);

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
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) listResults->Add(pDeltaOS[i][j]);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) listResults->Add(pDeltaSS[i][j]);
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
  UInt_t mask = handler->IsEventSelected();
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
  // Pair
  //----------------------------
  if (!PairTrkTrk()) return;
  hEvtCount->Fill("Pair Tracks",1);
  if (fDebug) Printf("Pair V0 & Trk done!");
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
    double phi  = track->Phi();
    double  pt  = track->Pt();
    double eta  = track->Eta();
    int charge  = track->Charge();
    int     id  = track->GetID();
    int  nhits  = track->GetTPCNcls();
    double dedx = track->GetTPCsignal();
    
    hPt->Fill(pt);
    hEta->Fill(eta);
    hPhi->Fill(phi);
    hNhits->Fill(nhits);
    h2PDedx->Fill(track->P()*charge, dedx);

    //DCA Cut
    double dcaxy = -999, dcaz = -999;
    if(!GetDCA(dcaxy,dcaz,track)) continue;
    if (fabs(dcaxy) > fDcaCutXY) continue;
    if (fabs(dcaz) > fDcaCutZ) continue;
    hDcaXY->Fill(fabs(dcaxy));
    hDcaZ->Fill(fabs(dcaz));

    int pid = GetPIDofParticle(track);
    if (isnan(pid)) continue;
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

    int label = track->GetLabel();
    vecParticle.emplace_back(std::array<double,5>{pt,eta,phi,(double)pid,(double)label});
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::PairTrkTrk()
{
  for (auto particle_1 : vecParticle) {
    double pt_1 = particle_1[0];
    double eta_1 = particle_1[1];
    double phi_1 = particle_1[2];
    int id_1 = particle_1[3];
    int label_1 = particle_1[4];

    for (auto particle_2 : vecParticle) {
      double pt_2 = particle_2[0];
      double eta_2 = particle_2[1];
      double phi_2 = particle_2[2];
      int id_2 = particle_2[3];
      int label_2 = particle_2[4];

      if(label_1 == label_2) continue;
      double delta = cos(phi_1 - phi_2);

      //id    211  321 2212
      //211    00   01   02
      //321    10   11   12
      //2212   20   21   22

      int pid_row = 999;
      int pid_col = 999;
      if(abs(id_1) == 211 && abs(id_2) == 211)        {pid_row = 0; pid_col = 0;}
      else if(abs(id_1) ==  211 && abs(id_2) ==  321) {pid_row = 0; pid_col = 1;}
      else if(abs(id_1) ==  211 && abs(id_2) == 2212) {pid_row = 0; pid_col = 2;}
      else if(abs(id_1) ==  321 && abs(id_2) ==  211) {pid_row = 1; pid_col = 0;}
      else if(abs(id_1) ==  321 && abs(id_2) ==  321) {pid_row = 1; pid_col = 1;}
      else if(abs(id_1) ==  321 && abs(id_2) == 2212) {pid_row = 1; pid_col = 2;}
      else if(abs(id_1) == 2212 && abs(id_2) ==  211) {pid_row = 2; pid_col = 0;}
      else if(abs(id_1) == 2212 && abs(id_2) ==  321) {pid_row = 2; pid_col = 1;}
      else if(abs(id_1) == 2212 && abs(id_2) == 2212) {pid_row = 2; pid_col = 2;}
      else continue;

      if(id_1 * id_2 < 0) pDeltaOS[pid_row][pid_col]->Fill(fCent,delta);
      else if(id_1 * id_2 > 0) pDeltaSS[pid_row][pid_col]->Fill(fCent,delta);
      else continue;
    }
  }
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskCMEPIDCVE::ResetVectors()
{
  vecParticle.clear();
  std::vector<std::array<double,5>>().swap(vecParticle);
}

//---------------------------------------------------

bool AliAnalysisTaskCMEPIDCVE::RejectEvtTFFit()
{
  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk) {
          delete aodTrk;
          continue;
      }
      if (aodTrk->TestFilterBit(32)) {
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  h2MultCent[0]->Fill(fCentV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t  multV0a = aodV0->GetMTotV0A();
  Float_t  multV0c = aodV0->GetMTotV0C();
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (fCentSPD0 < fCenCutLowPU->Eval(fCentV0M)) return false;
  if (fCentSPD0 > fCenCutHighPU->Eval(fCentV0M)) return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(fCentV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  if (isTightPileUp == true){
    Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
    Float_t nclsDif = Float_t(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot);
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