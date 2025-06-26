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
// CVE and PID CME analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

#include <TGraph.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TString.h>
#include <sys/time.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <array>
#include <unordered_map>
#include <vector>
// ROOT classes
#include "AliLog.h"
#include "TBits.h"
#include "TChain.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpline.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODVZERO.h"
#include "AliAODVertex.h"
// Alice classes
#include "AliEventCuts.h"
#include "AliEventplane.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliOADBContainer.h"
// Alice "V" classes
// Alice PID classes
#include "AliAnalysisTaskCVEPIDCMEDiff.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskCVEPIDCMEDiff);

namespace {
constexpr float LAMBDAMASS = 1.115683;
constexpr float MASSCUT    = 0.02;
constexpr int MASSBIN      = 60;
constexpr float INTGBIN    = 0.5;

static const std::array<TString, 125> runNumList18q{
  "296623", "296622", "296621", "296619", "296618", "296616", "296615", "296594", "296553", "296552", "296551", "296550", "296548", "296547",
  "296516", "296512", "296511", "296510", "296509", "296472", "296433", "296424", "296423", "296420", "296419", "296415", "296414", "296383",
  "296381", "296380", "296379", "296378", "296377", "296376", "296375", "296312", "296309", "296304", "296303", "296280", "296279", "296273",
  "296270", "296269", "296247", "296246", "296244", "296243", "296242", "296241", "296240", "296198", "296197", "296196", "296195", "296194",
  "296192", "296191", "296143", "296142", "296135", "296134", "296133", "296132", "296123", "296074", "296066", "296065", "296063", "296062",
  "296060", "296016", "295942", "295941", "295937", "295936", "295913", "295910", "295909", "295861", "295860", "295859", "295856", "295855",
  "295854", "295853", "295831", "295829", "295826", "295825", "295822", "295819", "295818", "295816", "295791", "295788", "295786", "295763",
  "295762", "295759", "295758", "295755", "295754", "295725", "295723", "295721", "295719", "295718", "295717", "295714", "295712", "295676",
  "295675", "295673", "295668", "295667", "295666", "295615", "295612", "295611", "295610", "295589", "295588", "295586", "295585"};
static const std::array<TString, 89> runNumList18r{
  "297595", "297590", "297588", "297558", "297544", "297542", "297541", "297540", "297537", "297512", "297483", "297479", "297452", "297451",
  "297450", "297446", "297442", "297441", "297415", "297414", "297413", "297406", "297405", "297380", "297379", "297372", "297367", "297366",
  "297363", "297336", "297335", "297333", "297332", "297317", "297311", "297310", "297278", "297222", "297221", "297218", "297196", "297195",
  "297193", "297133", "297132", "297129", "297128", "297124", "297123", "297119", "297118", "297117", "297085", "297035", "297031", "296966",
  "296941", "296938", "296935", "296934", "296932", "296931", "296930", "296903", "296900", "296899", "296894", "296852", "296851", "296850",
  "296848", "296839", "296838", "296836", "296835", "296799", "296794", "296793", "296790", "296787", "296786", "296785", "296784", "296781",
  "296752", "296694", "296693", "296691", "296690"};

std::array<double, 8> parV0{43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
std::array<double, 6> parV0CL0{0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
std::array<double, 8> parFB32{2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};

template <typename G>
inline float EvalGraph(G* g, float pt) {
  return (g ? g->Eval(pt) : 1.f);
}

static std::unordered_map<std::string, std::vector<double>> pt_ranges = {
    {"hadron", {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0}},
    {"proton", {0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0}},
    {"lambda", {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0}}
};

float GetNUACorThisPDG(TH2F* h2_pt_nua, float pt, float phi) {
  int bin = h2_pt_nua->FindBin(pt,phi);
  float nua = h2_pt_nua->GetBinContent(bin);
  return nua > 1e-6f ? nua : 0.f;
}

} // namespace

//---------------------------------------------------
AliAnalysisTaskCVEPIDCMEDiff::AliAnalysisTaskCVEPIDCMEDiff() : AliAnalysisTaskSE() {
}

//---------------------------------------------------
AliAnalysisTaskCVEPIDCMEDiff::AliAnalysisTaskCVEPIDCMEDiff(const char* name) : AliAnalysisTaskSE(name) {
  DefineInput(0, TChain::Class());
  DefineInput(1, TList::Class()); // NUE
  DefineInput(2, TList::Class()); // NUA
  DefineInput(3, TList::Class()); // TPC
  DefineInput(4, TList::Class()); // VZERO
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//------------------------------------------------

AliAnalysisTaskCVEPIDCMEDiff::~AliAnalysisTaskCVEPIDCMEDiff() {
  // Destructor
  // histograms are in the output list and deleted when the output
  if (runNumList) delete runNumList;
  runNumList = nullptr;
  if (fQAList) delete fQAList;
  fQAList = nullptr;
  if (fResultsList) delete fResultsList;
  fResultsList = nullptr;
  if (fListNUENUA) delete fListNUENUA;
  fListNUENUA = nullptr;
  if (fListPlaneNUA) delete fListPlaneNUA;
  fListPlaneNUA = nullptr;
  if (fListVZERO) delete fListVZERO;
  fListVZERO = nullptr;
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::Terminate(Option_t*) {
  // Terminate loop
  AliInfo("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::UserCreateOutputObjects() {
  /////////////////////////
  // Input stream
  /////////////////////////
  fListNUENUA     = dynamic_cast<TList*>(GetInputData(1));
  fListPlaneNUA   = dynamic_cast<TList*>(GetInputData(2));
  fListTPC        = dynamic_cast<TList*>(GetInputData(3));
  fListVZERO      = dynamic_cast<TList*>(GetInputData(4));

  /////////////////////////
  // Deal with calc flags
  /////////////////////////
  std::array<bool, 4> flags = {isCalculateLambdaHadron, isCalculateLambdaPion, isCalculateLambdaProton, isCalculateLambdaLambda};

  int enabledCount = std::count_if(flags.begin(), flags.end(), [](bool b) { return b; });

  if (enabledCount > 1) {
    AliFatal("Sorry but please do them one by one!");
  }

  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Rihan 18q/r Pile-up function
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fSPDCutPU = std::unique_ptr<TF1>(new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000));

    fV0CutPU = std::unique_ptr<TF1>(new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000));
    fV0CutPU->SetParameters(parV0.data());

    fCenCutLowPU = std::unique_ptr<TF1>(new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100));
    fCenCutLowPU->SetParameters(parV0CL0.data());

    fCenCutHighPU = std::unique_ptr<TF1>(new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100));
    fCenCutHighPU->SetParameters(parV0CL0.data());

    fMultCutPU = std::unique_ptr<TF1>(new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90));
    fMultCutPU->SetParameters(parFB32.data());
  } else {
    AliFatal("Sorry only support LHC18q/r dataset!");
  }

  ////////////////////////
  // NUE NUA
  ////////////////////////
  // Load Calibration Files
  // The global read-in lists and hists are loaded here.
  // They do not need to be loaded run by run.
  if (isDoNUE || isDoNUA || isDoLambdaNUE || isDoLambdaNUA) {
    if (!fListNUENUA) {
      AliFatal("NUE list not found");
      return;
    }
  }

  ////////////////////////
  /// TPC
  ////////////////////////
  if (fPlaneEstimator.EqualTo("TPC")) {
    if (!fListPlaneNUA || !fListNUENUA) {
        AliFatal("TPC NUA or NUE list not found");
        return;
    }
    if (isRecentreTPC) {
      if (!fListTPC) {
          AliFatal("TPC recentre list not found");
          return;
      }
      fQxTPCRunCentVz = (TProfile3D*)fListTPC->FindObject("QxTPCRunCentVz");
      fQyTPCRunCentVz = (TProfile3D*)fListTPC->FindObject("QyTPCRunCentVz");
      if (!fQxTPCRunCentVz || !fQyTPCRunCentVz) {
          AliFatal("TPC recentre hist not found");
          return;
      }
    }
  }

  ////////////////////////
  // VZERO
  ////////////////////////
  if (fPlaneEstimator.EqualTo("V0C")) {
    if (!fListVZERO) {
      AliFatal("V0C calibration list not found");
    }
    contQxncm = (AliOADBContainer*)fListVZERO->FindObject(Form("fqxc%im", 2)); // V0C Qx Mean
    contQyncm = (AliOADBContainer*)fListVZERO->FindObject(Form("fqyc%im", 2)); // V0C Qy Mean
  }

  //------------------
  // QA
  //------------------
  fQAList = new TList();
  fQAList->SetName("fQAList");
  fQAList->SetOwner(kTRUE);
  // event-wise
  fEvtCount = new TH1D("EvtCount", "Event Count", 20, 1, 21);
  fEvtCount->GetXaxis()->SetBinLabel(1, "All");
  fEvtCount->GetXaxis()->SetBinLabel(2, "Manager");
  fEvtCount->GetXaxis()->SetBinLabel(3, "Handler");
  fEvtCount->GetXaxis()->SetBinLabel(4, "fAOD");
  fEvtCount->GetXaxis()->SetBinLabel(5, "fPID");
  fEvtCount->GetXaxis()->SetBinLabel(6, "fMultSel");
  fEvtCount->GetXaxis()->SetBinLabel(7, "Trigger");
  fEvtCount->GetXaxis()->SetBinLabel(8, "Run Number");
  fEvtCount->GetXaxis()->SetBinLabel(9, "Calib Readin");
  fEvtCount->GetXaxis()->SetBinLabel(10, "Vertex");
  fEvtCount->GetXaxis()->SetBinLabel(11, "Centrality");
  fEvtCount->GetXaxis()->SetBinLabel(12, "Pile up");
  fEvtCount->GetXaxis()->SetBinLabel(13, "Reset Vector");
  fEvtCount->GetXaxis()->SetBinLabel(14, "Loop Track");
  fEvtCount->GetXaxis()->SetBinLabel(15, "Get Plane");
  fEvtCount->GetXaxis()->SetBinLabel(16, "Resolution");
  fEvtCount->GetXaxis()->SetBinLabel(17, "Loop V0");
  fEvtCount->GetXaxis()->SetBinLabel(18, "Pair V0Trk");
  fEvtCount->GetXaxis()->SetBinLabel(19, "Pair V0V0");
  fEvtCount->GetXaxis()->SetBinLabel(20, "Final Event");

  fQAList->Add(fEvtCount);

  ////////////////////////
  // Run Number Info
  ////////////////////////
  runNumList = new std::map<int, int>;

  if (fPeriod.EqualTo("LHC18q"))
    for (int i = 0; i < 125; i++) runNumList->insert(std::pair<int, int>(runNumList18q[i].Atoi(), i + 1));
  else if (fPeriod.EqualTo("LHC18r"))
    for (int i = 0; i < 89; i++) runNumList->insert(std::pair<int, int>(runNumList18r[i].Atoi(), i + 1));
  else
    return;
  fHistRunNumBin = new TH1I("runNumBin", "", (int)runNumList->size(), 1, (int)runNumList->size() + 1);
  std::map<int, int>::iterator iter;
  for (auto runNum : *runNumList) fHistRunNumBin->GetXaxis()->SetBinLabel(runNum.second, Form("%i", runNum.first));
  fQAList->Add(fHistRunNumBin);

  // Event-wise QA
  fHistCent[0]        = new TH1D("fHistCentBfCut", "Dist. of Centrality Before Cut", 80, 0., 80.);
  fHistCent[1]        = new TH1D("fHistCentAfCut", "Dist. of Centrality After Cut", 80, 0., 80.);
  fHistVz[0]          = new TH1D("fHistVzBfCut", "Dist of Vz Before Cut", 200, -20., 20.);
  fHistVz[1]          = new TH1D("fHistVzAfCut", "Dist of Vz After Cut", 200, -20., 20.);
  fHist2CentQA[0]     = new TH2D("fHist2CentQA_V0M_SPD1_BfCut", ";centV0M;centSPD1", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[1]     = new TH2D("fHist2CentQA_V0M_SPD1_AfCut", ";centV0M;centSPD1", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[2]     = new TH2D("fHist2CentQA_V0M_SPD0_BfCut", ";centV0M;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[3]     = new TH2D("fHist2CentQA_V0M_SPD0_AfCut", ";centV0M;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[4]     = new TH2D("fHist2CentQA_SPD1_SPD0_BfCut", ";centSPD1;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[5]     = new TH2D("fHist2CentQA_SPD1_SPD0_AfCut", ";centSPD1;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2MultCentQA[0] = new TH2D("fHist2MultCentQA_BfCut", ";centV0M;multFB32", 80, 0, 80., 50, 0, 5000.);
  fHist2MultCentQA[1] = new TH2D("fHist2MultCentQA_AfCut", ";centV0M;multFB32", 80, 0, 80., 50, 0, 5000.);

  for (auto& h : fHistCent) fQAList->Add(h);
  for (auto& h : fHistVz) fQAList->Add(h);
  for (auto& h : fHist2CentQA) fQAList->Add(h);
  for (auto& h : fHist2MultCentQA) fQAList->Add(h);

  // track-wise QA
  fHistPt     = new TH1D("fHistPt", ";p_{T}", 200, 0., 20.);
  fHistEta    = new TH1D("fHistEta", ";#eta", 100, -2., 2.);
  fHistNhits  = new TH1D("fHistNhits", ";nhits", 200, 0., 200.);
  fHist2PDedx = new TH2D("fHist2PDedx", ";pdedx", 400, -10., 10., 400, 0, 1000);
  fHistDcaXY  = new TH1D("fHistDcaXY", ";DcaXY", 500, 0., 1);
  fHistDcaZ   = new TH1D("fHistDcaZ", ";DcaZ", 500, 0., 1);
  fHistPhi[0] = new TH2D("fHistPhi_bfNUA", "hadron phi before NUA;p_{T};#phi", pt_ranges["hadron"].size() - 1, pt_ranges["hadron"].data(), 100, 0, TMath::TwoPi());
  fHistPhi[1] = new TH2D("fHistPhi_afNUA", "hadron phi after NUA;p_{T};#phi", pt_ranges["hadron"].size() - 1, pt_ranges["hadron"].data(), 100, 0, TMath::TwoPi());
  fQAList->Add(fHistPt);
  fQAList->Add(fHistEta);
  fQAList->Add(fHistNhits);
  fQAList->Add(fHist2PDedx);
  fQAList->Add(fHistDcaXY);
  fQAList->Add(fHistDcaZ);
  for (auto& h : fHistPhi) fQAList->Add(h);

  // Proton QA
  fHistProtonPt      = new TH1D("fHistProtonPt", "fHistProtonPt;p_{T}", 200, 0., 5.);
  fHistProtonEta     = new TH1D("fHistProtonEta", "fHistProtonEta;#eta", 100, -2., 2.);
  fHistProtonPhi[0]    = new TH2D("fHistProtonPhi_bfNUA", "Proton phi before NUA;p_{T};#phi", pt_ranges["proton"].size() - 1, pt_ranges["proton"].data(), 180, 0, TMath::TwoPi());
  fHistProtonPhi[1]    = new TH2D("fHistProtonPhi_afNUA", "Proton phi after NUA;p_{T};#phi", pt_ranges["proton"].size() - 1, pt_ranges["proton"].data(), 180, 0, TMath::TwoPi());
  fHistProtonPtDcaXY = new TH2D("fHistProtonPtDcaXY", "fHistProtonPtDcaXY;Pt;DcaXY", 23, 0.4, 5.0, 300, 0., 2);
  fHistProtonPtDcaZ  = new TH2D("fHistProtonPtDcaZ", "fHistProtonPtDcaZ;Pt;DcaZ", 23, 0.4, 5.0, 300, 0., 2);
  fQAList->Add(fHistProtonPt);
  fQAList->Add(fHistProtonEta);
  fQAList->Add(fHistProtonPhi[0]);
  fQAList->Add(fHistProtonPhi[1]);
  fQAList->Add(fHistProtonPtDcaXY);
  fQAList->Add(fHistProtonPtDcaZ);

  fHistAntiProtonPt      = new TH1D("fHistAntiProtonPt", "fHistAntiProtonPt;p_{T}", 200, 0., 5.);
  fHistAntiProtonEta     = new TH1D("fHistAntiProtonEta", "fHistAntiProtonEta;#eta", 100, -2., 2.);
  fHistAntiProtonPhi[0]  = new TH2D("fHistAntiProtonPhi_bfNUA", "AntiProton phi before NUA; p_{T}; #phi", pt_ranges["proton"].size() - 1, pt_ranges["proton"].data(), 180, 0, TMath::TwoPi());
  fHistAntiProtonPhi[1]  = new TH2D("fHistAntiProtonPhi_afNUA", "AntiProton phi after NUA; p_{T}; #phi", pt_ranges["proton"].size() - 1, pt_ranges["proton"].data(), 180, 0, TMath::TwoPi());
  fHistAntiProtonPtDcaXY = new TH2D("fHistAntiProtonPtDcaXY", "fHistAntiProtonPtDcaXY;Pt;DcaXY", 23, 0.4, 5.0, 300, 0., 2);
  fHistAntiProtonPtDcaZ  = new TH2D("fHistAntiProtonPtDcaZ", "fHistAntiProtonPtDcaZ;Pt;DcaZ", 23, 0.4, 5.0, 300, 0., 2);
  fQAList->Add(fHistAntiProtonPt);
  fQAList->Add(fHistAntiProtonEta);
  fQAList->Add(fHistAntiProtonPhi[0]);
  fQAList->Add(fHistAntiProtonPhi[1]);
  fQAList->Add(fHistAntiProtonPtDcaXY);
  fQAList->Add(fHistAntiProtonPtDcaZ);

  // V0s QA
  fHistV0Pt              = new TH1D("hV0Pt", "", 200, 0., 20.);
  fHistV0Eta             = new TH1D("hV0Eta", "", 200, -3., 3.);
  fHistV0DcatoPrimVertex = new TH1D("hV0DcaToPrimVertex", "", 200, 0., 20.);
  fHistV0CPA             = new TH1D("hV0CPA", "", 1000, 0.9, 1.);
  fHistV0DecayLength     = new TH1D("hV0DecayLength", "", 500, 0, 500.);
  fHistV0NegDaughterDca  = new TH1D("hV0NegDaughterDca", "", 200, 0., 20.);
  fHistV0PosDaughterDca  = new TH1D("hV0PosDaughterDca", "", 200, 0., 20.);
  fQAList->Add(fHistV0Pt);
  fQAList->Add(fHistV0Eta);
  fQAList->Add(fHistV0DcatoPrimVertex);
  fQAList->Add(fHistV0CPA);
  fQAList->Add(fHistV0DecayLength);
  fQAList->Add(fHistV0NegDaughterDca);
  fQAList->Add(fHistV0PosDaughterDca);

  // Lambda QA
  /// Lambda:
  fHistLambdaPt              = new TH1D("hLambdaPt", ";pT", 200, 0., 10.);
  fHistLambdaEta             = new TH1D("hLambdaEta", ";#eta", 200, -3., 3.);
  fHistLambdaPhi[0]           = new TH2D("fHistLambdaPhi_bfNUA","Lambda phi before NUA;p_{T};phi", pt_ranges["lambda"].size()-1, pt_ranges["lambda"].data(), 180, 0, TMath::TwoPi());
  fHistLambdaPhi[1]           = new TH2D("fHistLambdaPhi_afNUA","Lambda phi after NUA;p_{T};phi" , pt_ranges["lambda"].size()-1, pt_ranges["lambda"].data(), 180, 0, TMath::TwoPi());
  fHistLambdaDcaToPrimVertex = new TH1D("hLambdaDcaToPrimVertex", "DcatoPV", 200, 0., 20.);
  fHistLambdaNegDaughterDca  = new TH1D("hLambdaNegDaughterDca", "NegDaughterDcatoPV", 200, 0., 20.);
  fHistLambdaPosDaughterDca  = new TH1D("hLambdaPosDaughterDca", "PosDaughterDcatoPV", 200, 0., 20.);
  fHistLambdaCPA             = new TH1D("hLambdaCPA", ";cpa", 200, 0.9, 1.);
  fHistLambdaDecayLength     = new TH1D("hLambdaDecayLength", ";DecayLength", 250, 0., 500.);
  fHist2LambdaMassPtY        = new TH2D("fHist2LambdaMassPtY", ";pT;yH", 12, 0., 6., 20, -1., 1.);
  fQAList->Add(fHistLambdaPt);
  fQAList->Add(fHistLambdaEta);
  fQAList->Add(fHistLambdaPhi[0]);
  fQAList->Add(fHistLambdaPhi[1]);
  fQAList->Add(fHistLambdaDcaToPrimVertex);
  fQAList->Add(fHistLambdaNegDaughterDca);
  fQAList->Add(fHistLambdaPosDaughterDca);
  fQAList->Add(fHistLambdaCPA);
  fQAList->Add(fHistLambdaDecayLength);
  fQAList->Add(fHist2LambdaMassPtY);

  //AntiLambda
  fHistAntiLambdaPt              = new TH1D("hAntiLambdaPt", ";pT", 200, 0., 10.);
  fHistAntiLambdaEta             = new TH1D("hAntiLambdaEta", ";#eta", 200, -3., 3.);
  fHistAntiLambdaPhi[0]          = new TH2D("fHistAntiLambdaPhi_bfNUA","AntiLambda phi before NUA;p_{T};phi", pt_ranges["lambda"].size()-1, pt_ranges["lambda"].data(), 180, 0, TMath::TwoPi());
  fHistAntiLambdaPhi[1]          = new TH2D("fHistAntiLambdaPhi_afNUA","AntiLambda phi after NUA;p_{T};phi", pt_ranges["lambda"].size()-1, pt_ranges["lambda"].data(), 180, 0, TMath::TwoPi());
  fHistAntiLambdaDcaToPrimVertex = new TH1D("hAntiLambdaDcaToPrimVertex", "DcatoPV", 200, 0., 20.);
  fHistAntiLambdaNegDaughterDca  = new TH1D("hAntiLambdaNegDaughterDca", "NegDaughterDcatoPV", 200, 0., 20.);
  fHistAntiLambdaPosDaughterDca  = new TH1D("hAntiLambdaPosDaughterDca", "PosDaughterDcatoPV", 200, 0., 20.);
  fHistAntiLambdaCPA             = new TH1D("hAntiLambdaCPA", ";cpa", 200, 0.9, 1.);
  fHistAntiLambdaDecayLength     = new TH1D("hAntiLambdaDecayLength", ";DecayLength", 250, 0., 500.);
  fHist2AntiLambdaMassPtY        = new TH2D("fHist2AntiLambdaMassPtY", ";pT;yH", 12, 0., 6., 20, -1., 1.);
  fQAList->Add(fHistAntiLambdaPt);
  fQAList->Add(fHistAntiLambdaEta);
  fQAList->Add(fHistAntiLambdaPhi[0]);
  fQAList->Add(fHistAntiLambdaPhi[1]);
  fQAList->Add(fHistAntiLambdaDcaToPrimVertex);
  fQAList->Add(fHistAntiLambdaNegDaughterDca);
  fQAList->Add(fHistAntiLambdaPosDaughterDca);
  fQAList->Add(fHistAntiLambdaCPA);
  fQAList->Add(fHistAntiLambdaDecayLength);
  fQAList->Add(fHist2AntiLambdaMassPtY);

  fHist3LambdaCentPtMass      = new TH3D("fHist3LambdaCentPtMass", ";Centrality;Pt;Mass", 8, 0, 80, 20, 0, 10,
                                MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
  fQAList->Add(fHist3LambdaCentPtMass);

  fHist3LambdaCentPtMassWeighted = new TH3D("fHist3LambdaCentPtMassWeighted", ";Centrality;Pt;Mass", 8, 0, 80, 20, 0, 10,
                                MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
  fQAList->Add(fHist3LambdaCentPtMassWeighted);

  fHist3AntiLambdaCentPtMass = new TH3D("fHist3AntiLambdaCentPtMass", ";Centrality;Pt;Mass", 8, 0, 80, 20, 0, 10,
                                MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
  fQAList->Add(fHist3AntiLambdaCentPtMass);

  fHist3AntiLambdaCentPtMassWeighted = new TH3D("fHist3AntiLambdaCentPtMassWeighted", ";Centrality;Pt;Mass", 8, 0, 80, 20, 0, 10,
                                MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
  fQAList->Add(fHist3AntiLambdaCentPtMassWeighted);

  // Plane
  fProfile2DQxCentVz[0] = new TProfile2D("fProfile2DQxCentVz_bfRC", ";Centrality;VzBin;Qx", 70, 0, 70, 3, 0, 3);
  fProfile2DQyCentVz[0] = new TProfile2D("fProfile2DQyCentVz_bfRC", ";Centrality;VzBin;Qy", 70, 0, 70, 3, 0, 3);
  fQAList->Add(fProfile2DQxCentVz[0]);
  fQAList->Add(fProfile2DQyCentVz[0]);

  fProfile2DQxCentVz[1] = new TProfile2D("fProfile2DQxCentVz_afRC", ";Centrality;VzBin;Qx", 70, 0, 70, 3, 0, 3);
  fProfile2DQyCentVz[1] = new TProfile2D("fProfile2DQyCentVz_afRC", ";Centrality;VzBin;Qy", 70, 0, 70, 3, 0, 3);
  fQAList->Add(fProfile2DQxCentVz[1]);
  fQAList->Add(fProfile2DQyCentVz[1]);

  PostData(1, fQAList);
  if (fDebug) AliInfo("Post fQAList Data Success!");

  ////////////////////////
  // Results
  ////////////////////////
  fResultsList = new TList();
  fResultsList->SetName("fResultsList");
  fResultsList->SetOwner(kTRUE);
  fHist2Psi2 = new TH2D("fHist2Psi2", "Psi2 Dist.;centrality;#Psi2", 8, 0., 80., 100, 0., TMath::Pi());
  fResultsList->Add(fHist2Psi2);

  for (int iType = 0; iType < static_cast<int>(PairType::kNumPairs); iType++) {
    if (isCalculateLambdaHadron) {
      // Lambda - Hadron
      fHist3LambdaHadronMassIntg[iType] = new TH3D(Form("fHist3LambdaHadronMassIntg_%i", iType), Form("fHist3LambdaHadronMassIntg_%i", iType), 7, 0,
                                                   70, 1, 0, 1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fHist3LambdaHadronMassIntg[iType]);
      fProfile3DDeltaLambdaHadronMassIntg[iType] =
        new TProfile3D(Form("fProfile3DDeltaLambdaHadronMassIntg_%i", iType), Form("fProfile3DDeltaLambdaHadronMassIntg_%i", iType), 7, 0, 70, 1, 0,
                       1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DDeltaLambdaHadronMassIntg[iType]);
      fProfile3DGammaLambdaHadronMassIntg[iType] =
        new TProfile3D(Form("fProfile3DGammaLambdaHadronMassIntg_%i", iType), Form("fProfile3DGammaLambdaHadronMassIntg_%i", iType), 7, 0, 70, 1, 0,
                       1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DGammaLambdaHadronMassIntg[iType]);
    }
    if (isCalculateLambdaPion) {
      // Lambda - Pion
      fHist3LambdaPionMassIntg[iType] = new TH3D(Form("fHist3LambdaPionMassIntg_%i", iType), Form("fHist3LambdaPionMassIntg_%i", iType), 7, 0, 70, 1,
                                                 0, 1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fHist3LambdaPionMassIntg[iType]);
      fProfile3DDeltaLambdaPionMassIntg[iType] =
        new TProfile3D(Form("fProfile3DDeltaLambdaPionMassIntg_%i", iType), Form("fProfile3DDeltaLambdaPionMassIntg_%i", iType), 7, 0, 70, 1, 0, 1,
                       MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DDeltaLambdaPionMassIntg[iType]);
      fProfile3DGammaLambdaPionMassIntg[iType] =
        new TProfile3D(Form("fProfile3DGammaLambdaPionMassIntg_%i", iType), Form("fProfile3DGammaLambdaPionMassIntg_%i", iType), 7, 0, 70, 1, 0, 1,
                       MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DGammaLambdaPionMassIntg[iType]);
    }
    if (isCalculateLambdaLambda) {
      fHist3LambdaLambdaMassMass[iType] =
        new TH3D(Form("fHist3LambdaLambdaMassMass_%i", iType), Form("fHist3LambdaLambdaMassMass_%i", iType), 7, 0, 70, MASSBIN,
                 LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fHist3LambdaLambdaMassMass[iType]);
      fProfile3DDeltaLambdaLambdaMassMass[iType] =
        new TProfile3D(Form("fProfile3DDeltaLambdaLambdaMassMass_%i", iType), Form("fProfile3DDeltaLambdaLambdaMassMass_%i", iType), 7, 0, 70,
                       MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DDeltaLambdaLambdaMassMass[iType]);
      fProfile3DGammaLambdaLambdaMassMass[iType] =
        new TProfile3D(Form("fProfile3DGammaLambdaLambdaMassMass_%i", iType), Form("fProfile3DGammaLambdaLambdaMassMass_%i", iType), 7, 0, 70,
                       MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DGammaLambdaLambdaMassMass[iType]);
    }
    if (isCalculateLambdaProton) {
      // Lambda - Proton
      // Inv Mass
      fHist3LambdaProtonMassIntg[iType] = new TH3D(Form("fHist2LambdaProtonMassIntg_%i", iType), Form("fHist2LambdaProtonMassIntg_%i", iType), 7, 0,
                                                   70, 1, 0, 1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fHist3LambdaProtonMassSPt[iType]  = new TH3D(Form("fHist3LambdaProtonMassSPt_%i", iType), Form("fHist3LambdaProtonMassSPt_%i", iType), 7, 0, 70,
                                                   3, 0, 3, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fHist3LambdaProtonMassDEta[iType] = new TH3D(Form("fHist3LambdaProtonMassDEta_%i", iType), Form("fHist3LambdaProtonMassDEta_%i", iType), 7, 0,
                                                   70, 2, 0, 2, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fHist3LambdaProtonMassIntg[iType]);
      fResultsList->Add(fHist3LambdaProtonMassSPt[iType]);
      fResultsList->Add(fHist3LambdaProtonMassDEta[iType]);
      // Diff δ
      fProfile3DDeltaLambdaProtonMassIntg[iType] =
        new TProfile3D(Form("fProfile3DDeltaLambdaProtonMassIntg_%i", iType), Form("fProfile3DDeltaLambdaProtonMassIntg_%i", iType), 7, 0, 70, 1, 0,
                       1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fProfile3DDeltaLambdaProtonMassSPt[iType] =
        new TProfile3D(Form("fProfile3DDeltaLambdaProtonMassSPt_%i", iType), Form("fProfile3DDeltaLambdaProtonMassSPt_%i", iType), 7, 0, 70, 3, 0,
                       3, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fProfile3DDeltaLambdaProtonMassDEta[iType] =
        new TProfile3D(Form("fProfile3DDeltaLambdaProtonMassDEta_%i", iType), Form("fProfile3DDeltaLambdaProtonMassDEta_%i", iType), 7, 0, 70, 2, 0,
                       2, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DDeltaLambdaProtonMassIntg[iType]);
      fResultsList->Add(fProfile3DDeltaLambdaProtonMassSPt[iType]);
      fResultsList->Add(fProfile3DDeltaLambdaProtonMassDEta[iType]);
      //  γ
      fProfile3DGammaLambdaProtonMassIntg[iType] =
        new TProfile3D(Form("fProfile3DGammaLambdaProtonMassIntg_%i", iType), Form("fProfile3DGammaLambdaProtonMassIntg_%i", iType), 7, 0, 70, 1, 0,
                       1, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fProfile3DGammaLambdaProtonMassSPt[iType] =
        new TProfile3D(Form("fProfile3DGammaLambdaProtonMassSPt_%i", iType), Form("fProfile3DGammaLambdaProtonMassSPt_%i", iType), 7, 0, 70, 3, 0,
                       3, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fProfile3DGammaLambdaProtonMassDEta[iType] =
        new TProfile3D(Form("fProfile3DGammaLambdaProtonMassDEta_%i", iType), Form("fProfile3DGammaLambdaProtonMassDEta_%i", iType), 7, 0, 70, 2, 0,
                       2, MASSBIN, LAMBDAMASS - MASSCUT, LAMBDAMASS + MASSCUT);
      fResultsList->Add(fProfile3DGammaLambdaProtonMassIntg[iType]);
      fResultsList->Add(fProfile3DGammaLambdaProtonMassSPt[iType]);
      fResultsList->Add(fProfile3DGammaLambdaProtonMassDEta[iType]);
    }
  }
  PostData(2, fResultsList);
  if (fDebug) AliInfo("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::UserExec(Option_t*) {
  if (fDebug) AliInfo("===============================We are in UserExec!!!===================================");
  fEvtCount->Fill(1);

  bool isNeedPairV0Trk = isCalculateLambdaProton || isCalculateLambdaPion || isCalculateLambdaHadron;
  bool isNeedPairV0V0  = isCalculateLambdaLambda;

  // ----------------------------
  // Handle
  // ----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else {
    fEvtCount->Fill(2);
  }

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else
    fEvtCount->Fill(3);

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else
    fEvtCount->Fill(4);

  fPIDResponse = handler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else
    fEvtCount->Fill(5);

  fMultSel = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if (!fMultSel) {
    AliError(Form("%s: Could not get AliMultSelection", GetName()));
  } else
    fEvtCount->Fill(6);

  if (!manager || !handler || !fAOD || !fPIDResponse || !fMultSel) return;
  if (fDebug) AliInfo("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  unsigned int mask   = handler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB"))
    isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7"))
    isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kSemiCentral"))
    isTrigselected = mask & (AliVEvent::kINT7 | AliVEvent::kSemiCentral);
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
    isTrigselected = mask & (AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral);

  if (isTrigselected == false) return;

  fEvtCount->Fill(7);
  if (fDebug) AliInfo("trigger done!");

  //----------------------------
  // Run Number
  //----------------------------
  if (!runNumList || runNumList->empty()) {
    AliError("runNumList not initialized!");
    return;
  }
  fRunNum = fAOD->GetRunNumber();
  if (!(runNumList->find(fRunNum) != runNumList->end())) {
    AliError(Form("Run number %d not in runNumList! No calib files", fRunNum));
    return;
  }
  fEvtCount->Fill(8);

  //----------------------------
  // Reset event-by-event Vectors
  //----------------------------
  ResetVectors();
  fEvtCount->Fill(9);

  //----------------------------
  // Run-by-run Calib load
  //----------------------------
  if (fRunNum != fOldRunNum) {
    // Load the run dependent calibration hist
    if (!LoadCalibHistForThisRun()) {
        AliError(Form("loading calib hist for this run failed %d", fRunNum));
        return;
    }
    fRunNumBin = runNumList->at(fRunNum);
    fOldRunNum = fRunNum;
    if (fRunNumBin < 0) return;
  }
  fHistRunNumBin->Fill(fRunNumBin);
  fEvtCount->Fill(10);
  if (fDebug) AliInfo("calib file load done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  if (!fVtx) return;
  fVtx->GetXYZ(fVertex.data());
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (!vtSPD) return;
  if (fabs(fVertex[0]) < 1e-6 || fabs(fVertex[1]) < 1e-6 || fabs(fVertex[2]) < 1e-6) return;
  fHistVz[0]->Fill(fVertex[2]);
  if (fabs(fVertex[2]) > fVzCut) return;
  if (fVtx->GetNContributors() < 2 || vtSPD->GetNContributors() < 1) return;
  fHistVz[1]->Fill(fVertex[2]);
  for (int i = 0; i < 20; ++i) {
    if (fVertex[2] > -10 + i * 1 && fVertex[2] < -10 + (i + 1) * 1) {
      fVzBin = i;
      break;
    }
  }
  if (fVzBin < 0) return;
  fEvtCount->Fill(11);
  if (fDebug) AliInfo("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  float centV0M  = fMultSel->GetMultiplicityPercentile("V0M");
  float centSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
  float centSPD1 = fMultSel->GetMultiplicityPercentile("CL1");

  // we use centV0M as the default centrality
  fCent = centV0M;
  fHist2CentQA[0]->Fill(centV0M, centSPD1);
  fHist2CentQA[2]->Fill(centV0M, centSPD0);
  fHist2CentQA[4]->Fill(centSPD1, centSPD0);
  if (fabs(fCent - centSPD1) > 7.5) return;
  fHist2CentQA[1]->Fill(centV0M, centSPD1);
  fHist2CentQA[3]->Fill(centV0M, centSPD0);
  fHist2CentQA[5]->Fill(centSPD1, centSPD0);
  if (fCent < 0.f || fCent >= 70.f) return;

  // PF-Preview comment to make the shape edge of cent dist
  if (fCent > 30 && fCent < 50) {
    if (!(mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral))) return;
  } else {
    if (!(mask & AliVEvent::kINT7)) return;
  }

  fCentBin = (int)fCent / 10;
  fHistCent[0]->Fill(fCent);
  fEvtCount->Fill(12);
  // load cent dependent nue graph
  if (fCentBin != fOldCentBin) {
    if (!LoadNUENUAGraphForThisCent()) {
        AliError("NUE or NUA list loading failed");
        return;
    }
    fOldCentBin = fCentBin;
  }
  if (fDebug) AliInfo("centrality done!");

  //----------------------------
  // TPC recenter
  //----------------------------
  if (isRecentreTPC) {
      if (!fQxTPCRunCentVz || !fQyTPCRunCentVz) {
          fQxMeanTPC = 0.f;
          fQyMeanTPC = 0.f;
      } else {
          float vzBin = (fVertex[2] < -2.) ? 0.5 : (fVertex[2] > 2.) ? 2.5 : 1.5;
          int bin = fQxTPCRunCentVz->FindBin(fRunNumBin - 0.5, fCent, vzBin);
          if (bin < 1) {
              AliError("TPC QVector Calibration histogram bin out of range");
              fQxMeanTPC = 0.f;
              fQyMeanTPC = 0.f;
          } else {
              fQxMeanTPC = fQxTPCRunCentVz->GetBinContent(bin);
              fQyMeanTPC = fQyTPCRunCentVz->GetBinContent(bin);
          }
      }
  }
  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    if (!RejectEvtTFFit(centSPD0)) return;
  }
  fHistCent[1]->Fill(fCent);
  fEvtCount->Fill(13);
  if (fDebug) AliInfo("pile-up done!");
  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  // must loop tracks
  if (!LoopTracks()) return;
  fEvtCount->Fill(14);
  if (fDebug) AliInfo("Loop Tracks done!");
  //----------------------------
  // Plane
  //----------------------------
  if (fPlaneEstimator.EqualTo("TPC")) {
    fPsi2 = GetTPCPlane();
  } else if (fPlaneEstimator.EqualTo("V0C")) {
    fPsi2 = GetV0CPlane(centSPD1);
  } else {
    std::cout << "Wrong fPlaneEstimator!" << std::endl;
    return;
  }
  if (std::isnan(fPsi2) || std::isinf(fPsi2)) {
    AliError("fPsi2 of this event is NaN or Inf");
    return;
  }
  float vzBin = (fVertex[2] < -2.) ? 0.5 : (fVertex[2] > 2.) ? 2.5 : 1.5;
  if(fWgtMult > 1.e-6) {
    fProfile2DQxCentVz[0]->Fill(fCent, vzBin, fSumQ2x/fWgtMult);
    fProfile2DQyCentVz[0]->Fill(fCent, vzBin, fSumQ2y/fWgtMult);
    fProfile2DQxCentVz[1]->Fill(fCent, vzBin, fSumQ2x/fWgtMult - fQxMeanTPC);
    fProfile2DQyCentVz[1]->Fill(fCent, vzBin, fSumQ2y/fWgtMult - fQyMeanTPC);

    fEvtCount->Fill(15);
  }
  if (fDebug) AliInfo("Get Plane done!");
  //----------------------------
  // Fill Resolution
  //----------------------------
  fHist2Psi2->Fill(fCent, fPsi2);
  fEvtCount->Fill(16);
  //----------------------------
  // Get Lambda Vector
  //----------------------------
  if (!LoopV0s()) return;
  fEvtCount->Fill(17);
  if (fDebug) AliInfo("Get Lambda Vector done!");
  //----------------------------
  // Pair V0 Trk
  //----------------------------
  if (isNeedPairV0Trk) {
    if (!PairV0Trk()) return;
    fEvtCount->Fill(18);
  }
  if (fDebug) AliInfo("Pair V0 & Trk done!");
  //----------------------------
  // Pair V0 V0
  //----------------------------
  if (isNeedPairV0V0) {
    if (!PairV0V0()) return;
    fEvtCount->Fill(19);
  }
  if (fDebug) AliInfo("Pair V0 & V0 done!");
  //------------------
  // Post output data.
  //------------------
  fEvtCount->Fill(20);
  PostData(1, fQAList);
  PostData(2, fResultsList);
  if (fDebug) AliInfo("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::LoopTracks() {
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
    // DCA Cut
    float dcaxy = -999, dcaz = -999;
    if (!GetDCA(dcaxy, dcaz, track)) continue;
    if (fabs(dcaxy) > 2.4) continue;
    if (fabs(dcaz) > 3.2) continue;
    //------------------
    // NUE & NUA
    //------------------
    float phi  = track->Phi();
    if (phi < 0) phi += 2 * TMath::Pi();
    float pt   = track->Pt();
    float eta  = track->Eta();
    int charge = track->Charge();
    int id     = track->GetID();
    int nhits  = track->GetTPCNcls();
    float dedx = track->GetTPCsignal();

    fHistPt->Fill(pt);
    fHistEta->Fill(eta);
    fHistNhits->Fill(nhits);
    fHist2PDedx->Fill(track->P() * charge, dedx);
    fHistPhi[0]->Fill(pt, phi);

    if (fPlaneEstimator.EqualTo("TPC")) {
      if (pt > 0.2 && pt < 2.0) {
        // weight just for plane
        float weight = 1.;
        float wEff = GetPlaneNUECor(charge, pt);
        if (wEff < 0) continue;
        else weight *= wEff;

        float wAcc = GetPlaneNUACor(charge, phi, eta, fVertex[2]);
        if (wAcc < 0) continue;
        else weight *= wAcc;

        fHistPhi[1]->Fill(pt, phi, wAcc);
        // Do we need to set pT as weight for Better resolution?
        fSumQ2x += weight * TMath::Cos(2 * phi);
        fSumQ2y += weight * TMath::Sin(2 * phi);
        fWgtMult += weight;
        mapTPCTrksIDPhiWgt[id] = {phi,weight};
      }
    }

    if (fFilterBit == 768 && isNarrowDcaCuts768) {
      if (fabs(dcaz) > 2.0) continue;
      if (fabs(dcaxy) > 7.0 * (0.0026 + 0.005 / TMath::Power(pt, 1.01))) continue;
    }
    fHistDcaXY->Fill(fabs(dcaxy));
    fHistDcaZ->Fill(fabs(dcaz));

    bool isPIDTrkWeWant = false;
    if (isCalculateLambdaPion) {
      isPIDTrkWeWant = CheckPIDofParticle(track, 1); // 1=pion
      if (fSpecialHadronDCAzMax > 0) {
        isPIDTrkWeWant = isPIDTrkWeWant && (fabs(dcaz) < fSpecialHadronDCAzMax);
      }
      if (fPionMinPt > 0) {
        isPIDTrkWeWant = isPIDTrkWeWant && (pt > fPionMinPt);
      }
    }
    if (isCalculateLambdaProton) {
      isPIDTrkWeWant = CheckPIDofParticle(track, 3); // 3=proton
      isPIDTrkWeWant = isPIDTrkWeWant && (pt < 5.0 && pt > 0.4);
      if (isProtonCustomizedDCACut) {
        isPIDTrkWeWant = isPIDTrkWeWant && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035 / TMath::Power(pt, 1.1)));
      }
      if (fSpecialProtonDCAzMax > 0) {
        isPIDTrkWeWant = isPIDTrkWeWant && (fabs(dcaz) < fSpecialProtonDCAzMax);
      }
      if (fProtonMinPt > 0) {
        isPIDTrkWeWant = isPIDTrkWeWant && (pt > fProtonMinPt);
      }
    }
    if (isCalculateLambdaHadron) {
      isPIDTrkWeWant = CheckPIDofParticle(track, 0); // 0 = hadron
      if (fSpecialHadronDCAzMax > 0) {
        isPIDTrkWeWant = isPIDTrkWeWant && (fabs(dcaz) < fSpecialHadronDCAzMax);
      }
      if (fHadronMinPt > 0) {
        isPIDTrkWeWant = isPIDTrkWeWant && (pt > fHadronMinPt);
      }
    }
    if (!isPIDTrkWeWant) continue;

    int code         = 0;
    float pid_weight = 1.;
    if (charge > 0) {
      if (isCalculateLambdaProton)
        code = 2212;
      else if (isCalculateLambdaPion)
        code = 211;
      else
        code = 999;
      fHistProtonPt->Fill(pt);
      fHistProtonEta->Fill(eta);
      fHistProtonPhi[0]->Fill(pt, phi);
      fHistProtonPtDcaXY->Fill(pt, fabs(dcaxy));
      fHistProtonPtDcaZ->Fill(pt, fabs(dcaz));

      if (isDoNUA) {
        float nua_pid_weight = GetPIDNUACor(code, pt, phi);
        fHistProtonPhi[1]->Fill(pt, phi, nua_pid_weight);
        if (nua_pid_weight > 0) pid_weight *= nua_pid_weight;
      }
      if (isDoNUE) {
        float nue_pid_weight = GetPIDNUECor(code, pt);
        if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
      }

    } else {
      if (isCalculateLambdaProton)
        code = -2212;
      else if (isCalculateLambdaPion)
        code = -211;
      else
        code = -999;
      fHistAntiProtonPt->Fill(pt);
      fHistAntiProtonEta->Fill(eta);
      fHistAntiProtonPhi[0]->Fill(pt, phi);
      fHistAntiProtonPtDcaXY->Fill(pt, fabs(dcaxy));
      fHistAntiProtonPtDcaZ->Fill(pt, fabs(dcaz));

      if (isDoNUA) {
        float nua_pid_weight = GetPIDNUACor(code, pt, phi);
        fHistAntiProtonPhi[1]->Fill(pt, phi, nua_pid_weight);
        if (nua_pid_weight > 0) pid_weight *= nua_pid_weight;
      }
      if (isDoNUE) {
        float nue_pid_weight = GetPIDNUECor(code, pt);
        if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
      }
    }
    //vecParticle [pt,eta,phi,id,pdgcode,pidweight]
    vecParticle.emplace_back(std::array<float, 6>{pt, eta, phi, (float)id, (float)code, pid_weight});
  }

  if (fabs(fSumQ2x) < 1.e-6 || fabs(fSumQ2y) < 1.e-6 || fWgtMult < 1.e-5) return false;
  return true;
}

//---------------------------------------------------

float AliAnalysisTaskCVEPIDCMEDiff::GetTPCPlane() {
  double q2x = fSumQ2x/fWgtMult;
  double q2y = fSumQ2y/fWgtMult;
  //recentre
  if (isRecentreTPC) {
    q2x = q2x - fQxMeanTPC;
    q2y = q2y - fQyMeanTPC;
  }
  return GetEventPlane(q2x, q2y, 2);
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::LoopV0s() {
  int nV0s = fAOD->GetNumberOfV0s();
  for (int iV0 = 0; iV0 < nV0s; iV0++) {
    AliAODv0* v0 = fAOD->GetV0(iV0);
    if (!v0) {
      AliError(Form("%s: Could not get v0s", GetName()));
      continue;
    }
    // Basic kinematic variable
    float pt      = v0->Pt();
    float eta     = v0->PseudoRapV0();
    float dcaToPV = v0->DcaV0ToPrimVertex();              // DCA to Primary Vertex
    float cpa     = v0->CosPointingAngle(fVertex.data()); // cosine pointing angle
    float dl      = v0->DecayLengthV0(fVertex.data());
    float dcaNeg  = v0->DcaNegToPrimVertex();
    float dcaPos  = v0->DcaPosToPrimVertex();
    fHistV0Pt->Fill(pt);
    fHistV0Eta->Fill(eta);
    fHistV0DcatoPrimVertex->Fill(dcaToPV);
    fHistV0CPA->Fill(cpa);
    fHistV0DecayLength->Fill(dl);
    fHistV0NegDaughterDca->Fill(dcaNeg);
    fHistV0PosDaughterDca->Fill(dcaPos);
    int id = iV0; // seems like no ID/label for V0?
    // V0 cut
    if (!IsGoodV0(v0)) continue;
    // V0 daughters cut
    AliAODTrack* nTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    AliAODTrack* pTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    if (!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;
    float nDcaPV = v0->DcaNegToPrimVertex();
    float pDcaPV = v0->DcaPosToPrimVertex();
    if (nDcaPV < fDaughtersDCAToPrimVtxMin || pDcaPV < fDaughtersDCAToPrimVtxMin) continue;
    int code = GetLambdaCode(pTrack, nTrack);
    if (abs(code) != 3122) continue;

    float phi         = v0->Phi();
    if (phi < 0) phi += 2 * TMath::Pi();
    int id_daughter_1 = v0->GetPosID();
    int id_daughter_2 = v0->GetNegID();
    float rap         = v0->RapLambda();

    // if a daughter has been used as a daughter from Lambda before(It happends), we have to refuse a new one.
    // [pt,eta,phi,id,pdgcode,weight,mass,id1,id2]
    bool found = std::any_of(vecParticleV0.begin(), vecParticleV0.end(), [id_daughter_1, id_daughter_2](std::array<float, 9>& arr) {
      return ((int)arr[7] == id_daughter_1 || (int)arr[8] == id_daughter_2 || (int)arr[7] == id_daughter_2 || (int)arr[8] == id_daughter_1);
    });
    if (found) continue;

    float mass = (code == 3122) ? v0->MassLambda() : v0->MassAntiLambda();

    // lambda mass and mom cuts
    bool isInLambdaRange = true;
    isInLambdaRange      = isInLambdaRange && (pt > 0.5 && pt < 10);
    isInLambdaRange      = isInLambdaRange && (fabs(rap) < 0.5);
    isInLambdaRange      = isInLambdaRange && (fabs(mass - LAMBDAMASS) < MASSCUT);
    if (fLambdaMinPt > 0) {
      isInLambdaRange = isInLambdaRange && (pt > fLambdaMinPt);
    }
    if (!isInLambdaRange) continue;

    if (code == 3122) {
      fHistLambdaPt->Fill(pt);
      fHistLambdaEta->Fill(eta);
      fHistLambdaPhi[0]->Fill(pt, phi);
      fHistLambdaDcaToPrimVertex->Fill(dcaToPV);
      fHistLambdaNegDaughterDca->Fill(nDcaPV);
      fHistLambdaPosDaughterDca->Fill(pDcaPV);
      fHistLambdaCPA->Fill(cpa);
      fHistLambdaDecayLength->Fill(dl);
      fHist3LambdaCentPtMass->Fill(fCent, pt, mass);
      fHist2LambdaMassPtY->Fill(pt, rap);
    } else {
      fHistAntiLambdaPt->Fill(pt);
      fHistAntiLambdaEta->Fill(eta);
      fHistAntiLambdaPhi[0]->Fill(pt, phi);
      fHistAntiLambdaDcaToPrimVertex->Fill(dcaToPV);
      fHistAntiLambdaNegDaughterDca->Fill(nDcaPV);
      fHistAntiLambdaPosDaughterDca->Fill(pDcaPV);
      fHistAntiLambdaCPA->Fill(cpa);
      fHistAntiLambdaDecayLength->Fill(dl);
      fHist3AntiLambdaCentPtMass->Fill(fCent, pt, mass);
      fHist2AntiLambdaMassPtY->Fill(pt, rap);
    }

    // lambda NUA
    float weight = 1.;
    if (isDoLambdaNUA) {
      float nua_weight = GetPIDNUACor(code, pt, phi);
      // NUA QA
      code == 3122 ? fHistLambdaPhi[1]->Fill(pt, phi, nua_weight) : fHistAntiLambdaPhi[1]->Fill(pt, phi, nua_weight);
      if (nua_weight > 0.f) weight *= nua_weight;
    }

    // lambda NUE
    if (isDoLambdaNUE) {
      float nue_weight = GetPIDNUECor(code, pt);
      if (nue_weight > 0.f) weight *= nue_weight;
    }

    if (code == 3122) {
      fHist3LambdaCentPtMassWeighted->Fill(fCent, pt, mass, weight);
    } else {
      fHist3AntiLambdaCentPtMassWeighted->Fill(fCent, pt, mass, weight);
    }

    //vecParticleV0 [pt,eta,phi,id,pdgcode,pidweight, mass, id_daughter_1, id_daughter_2]
    vecParticleV0.emplace_back(std::array<float, 9>{pt, eta, phi, (float)id, (float)code, weight, mass, (float)id_daughter_1, (float)id_daughter_2});
  } // loop V0 end
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::PairV0Trk() {
  // Lambda - X
  for (auto lambda : vecParticleV0) {
    float pt_lambda  = lambda[0];
    float eta_lambda = lambda[1];
    float phi_lambda = lambda[2];
    // int    id_lambda     = (int)lambda[3];
    int code_lambda     = (int)lambda[4];
    float weight_lambda = lambda[5];
    float mass_lambda   = lambda[6];
    int id_daughter_1   = (int)lambda[7];
    int id_daughter_2   = (int)lambda[8];

    //vecParticle [pt,eta,phi,id,pdgcode,pidweight]
    for (auto particle : vecParticle) {
      float pt        = particle[0];
      float eta       = particle[1];
      float phi       = particle[2];
      int id          = (int)particle[3];
      int code        = (int)particle[4];
      float pidweight = particle[5];
      if (id == id_daughter_1 || id == id_daughter_2) continue;

      float psi2_forThisPair = std::numeric_limits<float>::quiet_NaN();;
      if (fPlaneEstimator.EqualTo("TPC")) {
        psi2_forThisPair = GetTPCPlaneNoAutoCorr({id_daughter_1, id_daughter_2, id});
      } else if (fPlaneEstimator.EqualTo("V0C")) {
        psi2_forThisPair = fPsi2;
      } else {
        AliError("fPlaneEstimator is not set properly.");
        return false;
      }

      if (std::isnan(psi2_forThisPair) || std::isinf(psi2_forThisPair)) continue;
      float delta = TMath::Cos(phi_lambda - phi);
      float gamma = TMath::Cos(phi_lambda + phi - 2 * psi2_forThisPair);

      int nBits = 4;
      if (isCalculateLambdaProton) {
        TBits bitsLambdaProtonPair(nBits);
        bitsLambdaProtonPair.SetBitNumber(0, code_lambda == 3122 && code == 2212);
        bitsLambdaProtonPair.SetBitNumber(1, code_lambda == 3122 && code == -2212);
        bitsLambdaProtonPair.SetBitNumber(2, code_lambda == -3122 && code == 2212);
        bitsLambdaProtonPair.SetBitNumber(3, code_lambda == -3122 && code == -2212);

        float sumPtBin    = GetSumPtBin(pt_lambda + pt);
        float deltaEtaBin = GetDeltaEtaBin(fabs(eta_lambda - eta));

        for (int iBits = 0; iBits < nBits; iBits++) {
          float weight_all = weight_lambda * pidweight;
          if (bitsLambdaProtonPair.TestBitNumber(iBits)) {
            fHist3LambdaProtonMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda);
            fProfile3DDeltaLambdaProtonMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda, delta, weight_all);
            fProfile3DGammaLambdaProtonMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda, gamma, weight_all);

            fHist3LambdaProtonMassSPt[iBits]->Fill(fCent, sumPtBin, mass_lambda);
            fProfile3DDeltaLambdaProtonMassSPt[iBits]->Fill(fCent, sumPtBin, mass_lambda, delta, weight_all);
            fProfile3DGammaLambdaProtonMassSPt[iBits]->Fill(fCent, sumPtBin, mass_lambda, gamma, weight_all);

            fHist3LambdaProtonMassDEta[iBits]->Fill(fCent, deltaEtaBin, mass_lambda);
            fProfile3DDeltaLambdaProtonMassDEta[iBits]->Fill(fCent, deltaEtaBin, mass_lambda, delta, weight_all);
            fProfile3DGammaLambdaProtonMassDEta[iBits]->Fill(fCent, deltaEtaBin, mass_lambda, gamma, weight_all);
          }
        }
      }

      if (isCalculateLambdaHadron) {
        TBits bitsLambdaHadronPair(nBits);
        bitsLambdaHadronPair.SetBitNumber(0, code_lambda == 3122 && code > 0);
        bitsLambdaHadronPair.SetBitNumber(1, code_lambda == 3122 && code < 0);
        bitsLambdaHadronPair.SetBitNumber(2, code_lambda == -3122 && code > 0);
        bitsLambdaHadronPair.SetBitNumber(3, code_lambda == -3122 && code < 0);
        for (int iBits = 0; iBits < nBits; iBits++) {
          float weight_all = weight_lambda * pidweight;
          if (bitsLambdaHadronPair.TestBitNumber(iBits)) {
            fHist3LambdaHadronMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda);
            fProfile3DDeltaLambdaHadronMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda, delta, weight_all);
            fProfile3DGammaLambdaHadronMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda, gamma, weight_all);
          }
        }
      }

      if (isCalculateLambdaPion) {
        TBits bitsLambdaPionPair(nBits);
        bitsLambdaPionPair.SetBitNumber(0, code_lambda == 3122 && code == 211);
        bitsLambdaPionPair.SetBitNumber(1, code_lambda == 3122 && code == -211);
        bitsLambdaPionPair.SetBitNumber(2, code_lambda == -3122 && code == 211);
        bitsLambdaPionPair.SetBitNumber(3, code_lambda == -3122 && code == -211);
        for (int iBits = 0; iBits < nBits; iBits++) {
          float weight_all = weight_lambda * pidweight;
          if (bitsLambdaPionPair.TestBitNumber(iBits)) {
            fHist3LambdaPionMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda);
            fProfile3DDeltaLambdaPionMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda, delta, weight_all);
            fProfile3DGammaLambdaPionMassIntg[iBits]->Fill(fCent, INTGBIN, mass_lambda, gamma, weight_all);
          }
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::PairV0V0() {
  const size_t nV0 = vecParticleV0.size();
  for (size_t i = 0; i < nV0; ++i) {
    const auto& lambda_a = vecParticleV0[i];

    double phi_a        = lambda_a[2];
    int id_a            = (int)lambda_a[3];
    int code_a          = (int)lambda_a[4];
    double weight_a     = lambda_a[5];
    double mass_a       = lambda_a[6];
    int id_a_daughter_1 = (int)lambda_a[7];
    int id_a_daughter_2 = (int)lambda_a[8];

    for (size_t j = 0; j < nV0; ++j) {
      const auto& lambda_b = vecParticleV0[j];

      double phi_b        = lambda_b[2];
      int id_b            = (int)lambda_b[3];
      int code_b          = (int)lambda_b[4];
      double weight_b     = lambda_b[5];
      double mass_b       = lambda_b[6];
      int id_b_daughter_1 = (int)lambda_b[7];
      int id_b_daughter_2 = (int)lambda_b[8];

      bool b_same_id = id_a == id_b;
      if (b_same_id) continue;

      bool b_shared_daughter = (id_a_daughter_1 == id_b_daughter_1);
      b_shared_daughter      = b_shared_daughter || (id_a_daughter_1 == id_b_daughter_2);
      b_shared_daughter      = b_shared_daughter || (id_a_daughter_2 == id_b_daughter_1);
      b_shared_daughter      = b_shared_daughter || (id_a_daughter_2 == id_b_daughter_2);
      if (b_shared_daughter) continue;

      float psi2_forThisPair = std::numeric_limits<float>::quiet_NaN();;
      if (fPlaneEstimator.EqualTo("TPC")) {
        psi2_forThisPair = GetTPCPlaneNoAutoCorr({id_a_daughter_1, id_a_daughter_2, id_b_daughter_1, id_b_daughter_2});
      } else if (fPlaneEstimator.EqualTo("V0C")) {
        psi2_forThisPair = fPsi2;
      } else {
        AliError("fPlaneEstimator is not set properly.");
        return false;
      }

      if (std::isnan(psi2_forThisPair) || std::isinf(psi2_forThisPair)) continue;

      int nBits = 4;
      TBits bitsLambdaLambdaPair(nBits);
      bitsLambdaLambdaPair.SetBitNumber(0, code_a == 3122 && code_b == 3122);
      bitsLambdaLambdaPair.SetBitNumber(1, code_a == 3122 && code_b == -3122);
      bitsLambdaLambdaPair.SetBitNumber(2, code_a == -3122 && code_b == 3122);
      bitsLambdaLambdaPair.SetBitNumber(3, code_a == -3122 && code_b == -3122);

      double weight_all = weight_a * weight_b;
      double delta      = TMath::Cos(phi_a - phi_b);
      double gamma      = TMath::Cos(phi_a + phi_b - 2 * psi2_forThisPair);

      for (int iBits = 0; iBits < nBits; iBits++) {
        if (bitsLambdaLambdaPair.TestBitNumber(iBits)) {
          fHist3LambdaLambdaMassMass[iBits]->Fill(fCent, mass_a, mass_b);
          fProfile3DDeltaLambdaLambdaMassMass[iBits]->Fill(fCent, mass_a, mass_b, delta, weight_all);
          fProfile3DGammaLambdaLambdaMassMass[iBits]->Fill(fCent, mass_a, mass_b, gamma, weight_all);
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::ResetVectors() {
  fSumQ2x  = 0.f, fSumQ2y  = 0.f;
  fWgtMult = 0.f;
  fQxMeanTPC = 0.f, fQyMeanTPC = 0.f;
  std::unordered_map<int, std::pair<float,float>>().swap(mapTPCTrksIDPhiWgt);
  std::vector<std::array<float, 6>>().swap(vecParticle);
  std::vector<std::array<float, 9>>().swap(vecParticleV0);
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::LoadCalibHistForThisRun() {
  // 18q/r NUA
  if (fPlaneEstimator.EqualTo("TPC")) {
    hCorrectNUAPos = (TH3F*)fListPlaneNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d", 0, fRunNum));
    hCorrectNUANeg = (TH3F*)fListPlaneNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d", 0, fRunNum));
    if (!hCorrectNUAPos) return false;
    if (!hCorrectNUANeg) return false;
  }
  // V0C
  if (fPlaneEstimator.EqualTo("V0C")) {
    // for recenter
    hQx2mV0C = ((TH1D*)contQxncm->GetObject(fRunNum));
    hQy2mV0C = ((TH1D*)contQyncm->GetObject(fRunNum));
    // for gain equalization
    fHCorrectV0ChWeghts = (TH2F*)fListVZERO->FindObject(Form("hWgtV0ChannelsvsVzRun%d", fRunNum));
    if (!hQx2mV0C) {
      AliError("Not found hQx2mV0C");
      return false;
    }
    if (!hQy2mV0C) {
      AliError("Not found hQy2mV0C");
      return false;
    }
    if (!fHCorrectV0ChWeghts) {
      AliError("Not found fHCorrectV0ChWeghts");
      return false;
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::RejectEvtTFFit(float centSPD0) {
  const int nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  const int nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  int nITSCls          = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  if (!aodTrkl) return false;

  const int nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const int nTracks = fAOD->GetNumberOfTracks();
  int multTrk       = 0;
  for (int it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
    if (!aodTrk) continue;
    if (aodTrk->TestFilterBit(32)) {
      if ((fabs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2)) multTrk++;
    }
  }

  fHist2MultCentQA[0]->Fill(fCent, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  if (!aodV0) return false;
  const float multV0a  = aodV0->GetMTotV0A();
  const float multV0c  = aodV0->GetMTotV0C();
  float multV0Tot      = multV0a + multV0c;
  const auto multV0aOn = aodV0->GetTriggerChargeA();
  const auto multV0cOn = aodV0->GetTriggerChargeC();
  auto multV0On        = multV0aOn + multV0cOn;

  // pile-up cuts with logging
  if (centSPD0 < fCenCutLowPU->Eval(fCent)) return false;
  if (centSPD0 > fCenCutHighPU->Eval(fCent)) return false;
  if (nITSCls > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (multTrk < fMultCutPU->Eval(fCent)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  if (isTightPileUp) {
    int tpcClsTot = fAOD->GetNumberOfTPCClusters();
    float nclsDif = tpcClsTot - (53182.6 + 113.326 * multV0Tot - 0.000831275 * multV0Tot * multV0Tot);
    if (nclsDif > 200000) // can be varied to 150000, 200000
      return false;
  }

  fHist2MultCentQA[1]->Fill(fCent, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::AcceptAODTrack(AliAODTrack* track) {
  //------------------
  // track cut
  //------------------
  float pt = track->Pt();
  if (pt < 0.2 || pt > 5.0) return false;
  float eta = track->Eta();
  if (fabs(eta) > 0.8) return false;
  int nhits = track->GetTPCNcls();
  if (nhits < fNclsCut) return false;
  float dedx = track->GetTPCsignal();
  if (dedx < 10.0) return false;
  float chi2 = track->Chi2perNDF();
  if (chi2 < fChi2Min) return false;
  if (chi2 > fChi2Max) return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck) {
  if (pidToCheck == 0) return kTRUE; //// Charge Particles do not need PID check

  if (!fPIDResponse) {
    AliInfo("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..

  float nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  float trkPtPID = ftrack->Pt();

  /// Pion =>
  if (pidToCheck == 1) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.5 && fabs(nSigTPC) <= fNSigmaTPC) return true;
    if (trkPtPID > 0.5 && fabs(nSigRMS) <= fNSigmaRMS) return true;
    return false;
  }
  /// Kaon =>
  else if (pidToCheck == 2) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.45 && fabs(nSigTPC) <= fNSigmaTPC) return true;
    if (trkPtPID > 0.45 && fabs(nSigRMS) <= fNSigmaRMS) return true;
    return false;
  }
  /// proton =>
  else if (pidToCheck == 3) { ///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    bool isProton = fabs(nSigRMS) < fNSigmaRMS;

    // if(isUsePionRejection) {
    //   float nSigTPCPion = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
    //   float nSigTOFPion = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    //   float nSigRMSPion = TMath::Sqrt(nSigTPCPion*nSigTPCPion + nSigTOFPion*nSigTOFPion);
    //   float mom = ftrack->P();
    //   bool isPassPionRejection = (mom > 0.5 && fabs(nSigRMSPion) > 3.);
    //   isPassPionRejection = isPassPionRejection || (mom < 0.5 && fabs(nSigTPCPion) > 3.);
    //   isProton = isProton && isPassPionRejection;
    // }
    return isProton;
  } else {
    AliInfo("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return false;
  }
}

//---------------------------------------------------
float AliAnalysisTaskCVEPIDCMEDiff::GetPlaneNUECor(int charge, float pt)
{
  if (charge == 0) return 1.f;
  if (!gNUEPosHadron_thisCent || !gNUENegHadron_thisCent) return 1.f;
  return (charge > 0) ? gNUEPosHadron_thisCent->Eval(pt) : gNUENegHadron_thisCent->Eval(pt);
}

//---------------------------------------------------

float AliAnalysisTaskCVEPIDCMEDiff::GetPlaneNUACor(int charge, float phi, float eta, float vz) {
  float weightNUA = 1;
  if (fVzBin < 0 || fCentBin < 0 || fRunNum < 0) return -1;
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) { // Rihan and Protty 's NUA Results
    if (charge > 0) {
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz, phi, eta);
      if (hCorrectNUAPos->GetBinContent(iBinNUA) > 0) weightNUA = (float)hCorrectNUAPos->GetBinContent(iBinNUA);
      return weightNUA;
    } else if (charge < 0) {
      if (!hCorrectNUANeg) return -1;
      int iBinNUA = hCorrectNUANeg->FindBin(vz, phi, eta);
      if (hCorrectNUANeg->GetBinContent(iBinNUA) > 0) weightNUA = (float)hCorrectNUANeg->GetBinContent(iBinNUA);
      return weightNUA;
    }
    // In Rihan and Protty 's NUA results, the phi distribution is independent on centrality and particle charge
  }
  return weightNUA;
}

//---------------------------------------------------
//
float AliAnalysisTaskCVEPIDCMEDiff::GetPIDNUECor(int pdg, float pt)
{
  if (!isDoNUE) return 1.f;

  switch (pdg) {
    case  999:  return EvalGraph(gNUEPosHadron_thisCent , pt);
    case -999:  return EvalGraph(gNUENegHadron_thisCent , pt);
    case  2212: return EvalGraph(gNUEProton_thisCent    , pt);
    case -2212: return EvalGraph(gNUEAntiProton_thisCent, pt);
    case  3122: return EvalGraph(gNUELambda_thisCent    , pt);
    case -3122: return EvalGraph(gNUEAntiLambda_thisCent, pt);
    default:    return 1.f;
  }
}

//---------------------------------------------------
//
float AliAnalysisTaskCVEPIDCMEDiff::GetPIDNUACor(int pdg, float pt, float phi)
{
  if (!isDoNUA) return 1.f;

  switch (pdg) {
    case  999:  return GetNUACorThisPDG(h2NUAPosHadron_thisCent , pt, phi);
    case -999:  return GetNUACorThisPDG(h2NUANegHadron_thisCent , pt, phi);
    case  2212: return GetNUACorThisPDG(h2NUAProton_thisCent    , pt, phi);
    case -2212: return GetNUACorThisPDG(h2NUAAntiProton_thisCent, pt, phi);
    case  3122: return GetNUACorThisPDG(h2NUALambda_thisCent    , pt, phi);
    case -3122: return GetNUACorThisPDG(h2NUAAntiLambda_thisCent, pt, phi);
    default:    return 1.f;
  }

}


//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::IsGoodV0(AliAODv0* aodV0) {
  // Offline reconstructed V0 only
  if (aodV0->GetOnFlyStatus()) return false;
  // Get daughters and check them
  AliAODTrack* myTrackPosTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
  AliAODTrack* myTrackNegTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
  if (!myTrackPosTest || !myTrackNegTest) {
    AliInfo("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return false;
  }
  // Unlike signs of daughters
  if (myTrackNegTest->Charge() * myTrackPosTest->Charge() > 0) return false;
  // Cosinus of pointing angle < 0.997
  float dCPA = aodV0->CosPointingAngle(fVertex.data());
  if (dCPA < fV0CPAMin) return false;
  // DCA of V0 < 1.5 cm
  float dV0Dca = aodV0->DcaV0ToPrimVertex();
  if (fabs(dV0Dca) > 1.5) return false;
  // V0 path length before decay 3-100 cm
  float dDecayLength = aodV0->DecayLengthV0(fVertex.data());
  if (dDecayLength > fV0DecayLengthMax) return false;
  if (dDecayLength < fV0DecayLengthMin) return false;
  // DCA between daughters < 0.5cm
  float dDCA = aodV0->DcaV0Daughters();
  if (dDCA > fV0DcaBetweenDaughtersMax) return false;
  return kTRUE;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::IsGoodDaughterTrack(const AliAODTrack* track) {
  // TPC refit
  if (!track->IsOn(AliAODTrack::kTPCrefit)) return false;
  // No kinks
  if (int(track->GetProdVertex()->GetType()) == AliAODVertex::kKink) return false;
  // Maximum value of transverse momentum
  float dPt = track->Pt();
  if (dPt > 20) return false;
  // Maximum value of pseudorapidity
  float dEta = track->Eta();
  if (fabs(dEta) > 0.8) return false;
  // Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2, 1);
  if (nCrossedRowsTPC < fDaughtersTPCNclsMin) return false;
  // Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return false;
  // [number of crossed rows] > 0.8 x [number of findable clusters].
  if (nCrossedRowsTPC / findable < fRatioCrossedRowsFindable) return false;
  return true;
}

//---------------------------------------------------

int AliAnalysisTaskCVEPIDCMEDiff::GetLambdaCode(const AliAODTrack* pTrack, const AliAODTrack* nTrack) {
  bool isLambda     = false;
  bool isAntiLambda = false;
  int code          = 0;

  // Λ-->(p+)+(π-)
  float nSigTPCPosProton = fabs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)); // TPC p+
  float nSigTPCNegPion   = fabs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));   // TPC π-
  //(Λ-)-->(p-)+(π+)
  float nSigTPCPosPion   = fabs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));   // TPC π+
  float nSigTPCNegProton = fabs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)); // TPC p-

  isLambda     = (nSigTPCPosProton < fDaughtersNSigmaTPC) && (nSigTPCNegPion < fDaughtersNSigmaTPC);
  isAntiLambda = (nSigTPCNegProton < fDaughtersNSigmaTPC) && (nSigTPCPosPion < fDaughtersNSigmaTPC);

  if (isLambda) code = 3122;
  if (isAntiLambda) code = -3122;
  if (isLambda && isAntiLambda) code = 0;
  return code;
}

//---------------------------------------------------

inline float AliAnalysisTaskCVEPIDCMEDiff::GetEventPlane(float qx, float qy, float harmonic) {
  float psi = (1. / harmonic) * TMath::ATan2(qy, qx);
  if (psi < 0) psi += TMath::TwoPi() / harmonic;
  return psi;
}

//---------------------------------------------------
float AliAnalysisTaskCVEPIDCMEDiff::GetTPCPlaneNoAutoCorr(std::vector<int> vec_id) {
  float psiNoAuto = std::numeric_limits<float>::quiet_NaN();;

  float tempSumQ2x  = fSumQ2x;
  float tempSumQ2y  = fSumQ2y;
  float tempWgtMult = fWgtMult;
  float repeQ2x = 0., repeQ2y = 0., repeWgtMult = 0.;

  std::vector<std::unordered_map<int, std::pair<float,float>>::iterator> vec_it;
  for (int id : vec_id) {
    vec_it.push_back(mapTPCTrksIDPhiWgt.find(id));
  }

  for (auto it : vec_it) {
    if (it != mapTPCTrksIDPhiWgt.end()) {
      float phi = it->second.first;
      float wgt = it->second.second;
      repeQ2x += wgt * TMath::Cos(2 * phi);
      repeQ2y += wgt * TMath::Sin(2 * phi);
      repeWgtMult += wgt;
    }
  }

  tempSumQ2x -= repeQ2x;
  tempSumQ2y -= repeQ2y;
  tempWgtMult -= repeWgtMult;

  if (tempWgtMult > 1.e-6) {
    //recentre
    tempSumQ2x /= tempWgtMult;
    tempSumQ2y /= tempWgtMult;
    if (isRecentreTPC) {
      tempSumQ2x -= fQxMeanTPC;
      tempSumQ2y -= fQyMeanTPC;
    }
    psiNoAuto = GetEventPlane(tempSumQ2x, tempSumQ2y, 2);
  } else {
    psiNoAuto = std::numeric_limits<float>::quiet_NaN();;
  }

  return psiNoAuto;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::GetDCA(float& dcaxy, float& dcaz, AliAODTrack* track) {
  if (!track) return false;
  double r[3];
  if (track->GetXYZ(r)) {
    dcaxy = r[0];
    dcaz  = r[1];
  } else {
    float dcax = r[0] - fVertex[0];
    float dcay = r[1] - fVertex[1];
    dcaz       = r[2] - fVertex[2];
    dcaxy      = sqrt(dcax * dcax + dcay * dcay);
  }
  return true;
}

//---------------------------------------------------
inline float AliAnalysisTaskCVEPIDCMEDiff::GetSumPtBin(float sumPt) {
  // 1~3 GeV return 0.5
  // 3~5 GeV return 1.5
  // 5~8 GeV return 2.5
  // else return -1
  if (sumPt > 1. && sumPt < 3.)
    return 0.5;
  else if (sumPt >= 3. && sumPt < 5.)
    return 1.5;
  else if (sumPt >= 5. && sumPt < 8.)
    return 2.5;
  else
    return -1.;
}

//---------------------------------------------------
inline float AliAnalysisTaskCVEPIDCMEDiff::GetDeltaEtaBin(float deltaEta) {
  //-0.6 ~ 0.6 return 0.5
  // < -0.6 & > 0.6 return 1.5
  // else return -1
  if (deltaEta < 0.6)
    return 0.5;
  else if (deltaEta > 0.6)
    return 1.5;
  else
    return -1.;
}

//---------------------------------------------------
inline float AliAnalysisTaskCVEPIDCMEDiff::GetDCABin(float dca) {
  dca = fabs(dca);
  // 0, 0.002 0.005 0.01
  if (dca > 0. && dca < 0.005)
    return 0.5;
  else if (dca > 0.005 && dca < 0.01)
    return 1.5;
  else if (dca > 0.01 && dca < 0.1)
    return 2.5;
  else if (dca > 0.1)
    return 3.5;
  else
    return -1.;
}

//---------------------------------------------------

float AliAnalysisTaskCVEPIDCMEDiff::GetV0CPlane(float centSPD1) {
  float psi2V0C = std::numeric_limits<float>::quiet_NaN();;
  float Qx2 = 0, Qy2 = 0, Mult = 0;

  // Loop Over VZERO Channels
  // Gain Equalization
  for (int iCh = 0; iCh < 32; ++iCh) { // just need C side
    // phi angle of the channel
    float phi = TMath::Pi() / 8. + TMath::Pi() / 4. * (iCh % 8);
    // energy of the channel
    float energy_raw   = 0.;
    AliAODVZERO* aodV0 = fAOD->GetVZEROData();
    energy_raw         = aodV0->GetMultiplicity(iCh);

    // get the gain equlization factor
    int ibinV0        = fHCorrectV0ChWeghts->FindBin(fVertex[2], iCh);
    float factor_gain = (float)fHCorrectV0ChWeghts->GetBinContent(ibinV0);

    float energy_gain = energy_raw * factor_gain;

    Qx2 += energy_gain * TMath::Cos(2 * phi);
    Qy2 += energy_gain * TMath::Sin(2 * phi);
    Mult += energy_gain;
  }
  if (Mult < 1.e-6) return std::numeric_limits<float>::quiet_NaN();;

  // VZERO Recenter
  //  get the mean Q vector from input file
  int iCentSPD  = (int)centSPD1;
  float Qx2Mean = hQx2mV0C->GetBinContent(iCentSPD + 1);
  float Qy2Mean = hQy2mV0C->GetBinContent(iCentSPD + 1);
  if (Qy2Mean < -900 || Qx2Mean < -900) return std::numeric_limits<float>::quiet_NaN();;

  // QA
  float vzBin = (fVertex[2] < -2.) ? 0.5 : (fVertex[2] > 2.) ? 2.5 : 1.5;
  fProfile2DQxCentVz[0]->Fill(fCent, vzBin, Qx2/Mult);
  fProfile2DQyCentVz[0]->Fill(fCent, vzBin, Qy2/Mult);

  // recenter the Q vector
  fSumQ2x = Qx2 - Qx2Mean;
  fSumQ2y = Qy2 - Qy2Mean;
  fWgtMult = Mult;

  fProfile2DQxCentVz[1]->Fill(fCent, vzBin, fSumQ2x/fWgtMult);
  fProfile2DQyCentVz[1]->Fill(fCent, vzBin, fSumQ2y/fWgtMult);

  psi2V0C = GetEventPlane(Qx2, Qy2, 2.);

  return psi2V0C;
}

//---------------------------------------------------
//
bool AliAnalysisTaskCVEPIDCMEDiff::LoadNUENUAGraphForThisCent() {
  TString period;
  if (fPeriod.EqualTo("LHC18q"))
    period = "18q";
  else if (fPeriod.EqualTo("LHC18r"))
    period = "18r";
  else
    return false;

  bool allFound = true;

  // --- Load NUE graphs ---
  if (isDoNUE) {
    std::vector<std::pair<TGraph**, TString>> nueGraphs = {
        {&gNUEPosHadron_thisCent, Form("nue_pt_poshadron_%s_cent%d", period.Data(), fCentBin)},
        {&gNUENegHadron_thisCent, Form("nue_pt_neghadron_%s_cent%d", period.Data(), fCentBin)},
        {&gNUEProton_thisCent,    Form("nue_pt_proton_%s_cent%d",    period.Data(), fCentBin)},
        {&gNUEAntiProton_thisCent,Form("nue_pt_antiproton_%s_cent%d",period.Data(), fCentBin)}
    };

    for (auto& p : nueGraphs) {
      *(p.first) = (TGraph*)fListNUENUA->FindObject(p.second);
      if (!*(p.first)) {
        AliError(Form("Could not find NUE graph: %s", p.second.Data()));
        allFound = false;
      }
    }
  }

  if (isDoLambdaNUE) {
    // Separately load Lambda NUE graphs
    gNUELambda_thisCent     = (TGraph*)fListNUENUA->FindObject(Form("nue_pt_lambda_%s_cent%d",     period.Data(), fCentBin));
    gNUEAntiLambda_thisCent = (TGraph*)fListNUENUA->FindObject(Form("nue_pt_antilambda_%s_cent%d", period.Data(), fCentBin));

    if (!gNUELambda_thisCent) {
      AliError("Could not find NUE graph: nue_pt_lambda");
      allFound = false;
    }
    if (!gNUEAntiLambda_thisCent) {
      AliError("Could not find NUE graph: nue_pt_antilambda");
      allFound = false;
    }
  }

  // --- Load NUA graphs ---
  if (isDoNUA) {
    std::vector<std::pair<TH2F**, TString>> nuaHist2 = {
        {&h2NUAPosHadron_thisCent, Form("nua_pt_poshadron_%s_cent%d", period.Data(), fCentBin)},
        {&h2NUANegHadron_thisCent, Form("nua_pt_neghadron_%s_cent%d", period.Data(), fCentBin)},
        {&h2NUAProton_thisCent,    Form("nua_pt_proton_%s_cent%d",    period.Data(), fCentBin)},
        {&h2NUAAntiProton_thisCent,Form("nua_pt_antiproton_%s_cent%d",period.Data(), fCentBin)}
    };

    for (auto& p : nuaHist2) {
      *(p.first) = (TH2F*)fListNUENUA->FindObject(p.second);
      if (!*(p.first)) {
        AliError(Form("Could not find NUA graph: %s", p.second.Data()));
        allFound = false;
      }
    }
  }

  if (isDoLambdaNUA) {
    // Separately load Lambda NUA graphs
    h2NUALambda_thisCent     = (TH2F*)fListNUENUA->FindObject(Form("nua_pt_lambda_%s_cent%d",     period.Data(), fCentBin));
    h2NUAAntiLambda_thisCent = (TH2F*)fListNUENUA->FindObject(Form("nua_pt_antilambda_%s_cent%d", period.Data(), fCentBin));

    if (!h2NUALambda_thisCent) {
      AliError("Could not find NUA graph: nua_pt_lambda");
      allFound = false;
    }
    if (!h2NUAAntiLambda_thisCent) {
      AliError("Could not find NUA graph: nua_pt_antilambda");
      allFound = false;
    }
  }

  return allFound;
}
