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
#include "AliAnalysisTaskCVEVZEROCalib.h"

ClassImp(AliAnalysisTaskCVEVZEROCalib);

//---------------------------------------------------
AliAnalysisTaskCVEVZEROCalib::AliAnalysisTaskCVEVZEROCalib() :
  AliAnalysisTaskSE(),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fAOD(nullptr),
  fUtils(nullptr),
  fMultSel(nullptr),
  fRunNum(-999),
  fOldRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCent(-999),
  fCentBin(-999),
  fCentV0M(-999),
  fCentTRK(-999),
  fCentSPD0(-999),
  fCentSPD1(-999),
  fVzBinNum(3),
  fCentBinNum(40),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fListVZEROCalib(nullptr),
  hMultV0Read(nullptr),
  hMultV0(nullptr),
  contMult(nullptr),
  contQxncm(nullptr),
  contQyncm(nullptr),
  contQxnam(nullptr),
  contQynam(nullptr),
  fHCorrectV0ChWeghts(nullptr),
  fQAList(nullptr),
  fEvtCount(nullptr),
  runNumList(nullptr),
  fHistRunNumBin(nullptr),
  fResultsList(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 3; i++) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; i++) pV0YMeanRead[i] = nullptr;
  for (int i = 0; i < 2; i++) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) hQy2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileChanalMult[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0ACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0ACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0CVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0CVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0AVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0AVz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0C_meanQ2x_cent_vz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0C_meanQ2y_cent_vz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0A_meanQ2x_cent_vz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0A_meanQ2y_cent_vz[i] = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskCVEVZEROCalib::AliAnalysisTaskCVEVZEROCalib(const char *name) :
  AliAnalysisTaskSE(name),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fAOD(nullptr),
  fUtils(nullptr),
  fMultSel(nullptr),
  fRunNum(-999),
  fOldRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCent(-999),
  fCentBin(-999),
  fCentV0M(-999),
  fCentTRK(-999),
  fCentSPD0(-999),
  fCentSPD1(-999),
  fVzBinNum(3),
  fCentBinNum(40),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fListVZEROCalib(nullptr),
  hMultV0Read(nullptr),
  hMultV0(nullptr),
  contMult(nullptr),
  contQxncm(nullptr),
  contQyncm(nullptr),
  contQxnam(nullptr),
  contQynam(nullptr),
  fHCorrectV0ChWeghts(nullptr),
  fQAList(nullptr),
  fEvtCount(nullptr),
  runNumList(nullptr),
  fHistRunNumBin(nullptr),
  fResultsList(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 3; i++) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; i++) pV0YMeanRead[i] = nullptr;
  for (int i = 0; i < 2; i++) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) hQy2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileChanalMult[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0ACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0ACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0CVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0CVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2xV0AVz[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileQ2yV0AVz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0C_meanQ2x_cent_vz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0C_meanQ2y_cent_vz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0A_meanQ2x_cent_vz[i] = nullptr;
  for (int i = 0; i < 150; i++) p2D_V0A_meanQ2y_cent_vz[i] = nullptr;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskCVEVZEROCalib::~AliAnalysisTaskCVEVZEROCalib()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fQAList) delete fQAList;
  if (fResultsList) delete fResultsList;
}

//---------------------------------------------------

void AliAnalysisTaskCVEVZEROCalib::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskCVEVZEROCalib::UserCreateOutputObjects()
{
  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Dobrin 15o pass2 Pile-up function
  if (fPeriod.EqualTo("LHC15o")) {
    fSPDCutPU = new TF1("fSPDCutPU", "450. + 3.9*x", 0, 50000);

    Double_t parV0[8] = {33.4237, 0.953516, 0.0712137, 227.923, 8.9239, -0.00319679, 0.000306314, -7.6627e-07};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    Double_t parV0CL0[6] = {0.0193587, 0.975914, 0.675714, 0.0292263, -0.000549509, 5.86421e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    Double_t parFB32[9] = {-812.822, 6.41796, 5421.83, -0.382601, 0.0299686, -26.6249, 321.388, -0.82615, 0.0167828};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
  }

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

  ////////////////////////
  // VZERO
  ////////////////////////
  if (!fListVZEROCalib) {
   std::cout<<("VZERO calibration list not found")<<std::endl;
   return;
  }
  if (fPeriod.EqualTo("LHC10h")) {
    // Read GE Hists (x_y) : (iCh_runnumBin)
    hMultV0Read = (TH2D*)fListVZEROCalib->FindObject("hMultV0");
    // Read Qx/y Mean Hists (x_y_z) : (runnumBin_centBin_vzBin)
    pV0XMeanRead[1] = (TProfile3D*)fListVZEROCalib->FindObject("pV0CCosMean");
    pV0YMeanRead[1] = (TProfile3D*)fListVZEROCalib->FindObject("pV0CSinMean");
    pV0XMeanRead[2] = (TProfile3D*)fListVZEROCalib->FindObject("pV0ACosMean");
    pV0YMeanRead[2] = (TProfile3D*)fListVZEROCalib->FindObject("pV0ASinMean");
  }
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    contQxncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxc%im",2)); // V0C Qx Mean
    contQyncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqyc%im",2)); // V0C Qy Mean
    contQxnam = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxa%im",2)); // V0A Qx Mean
    contQynam = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqya%im",2)); // V0A Qy Mean
    for (int i = 0; i < 2; i++) {
      hQx2mV0[i] = new TH1D();
      hQy2mV0[i] = new TH1D();
    }
  }
  //15 V0 Mult
  if (fPeriod.EqualTo("LHC15o")) {
    contMult = (AliOADBContainer*)fListVZEROCalib->FindObject("hMultV0BefCorPfpx");
    hMultV0  = new TH1D();
  }
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fHCorrectV0ChWeghts = new TH2F();
  }

  //------------------
  // QA
  //------------------
  fQAList = new TList();
  fQAList -> SetName("fQAList");
  fQAList -> SetOwner(kTRUE);
  // event-wise
  fEvtCount = new TH1D("EvtCount", "Event Count", 14, 1, 15);
  fEvtCount->GetXaxis()->SetBinLabel(1,"All");
  fEvtCount->GetXaxis()->SetBinLabel(2,"Manager");
  fEvtCount->GetXaxis()->SetBinLabel(3,"Handler");
  fEvtCount->GetXaxis()->SetBinLabel(4,"fAOD");
  fEvtCount->GetXaxis()->SetBinLabel(5,"fPID");
  fEvtCount->GetXaxis()->SetBinLabel(6,"fUtils");
  fEvtCount->GetXaxis()->SetBinLabel(7,"fMultSel");
  fEvtCount->GetXaxis()->SetBinLabel(8,"Trigger");
  fEvtCount->GetXaxis()->SetBinLabel(9,"Run Number");
  fEvtCount->GetXaxis()->SetBinLabel(10,"Read in");
  fEvtCount->GetXaxis()->SetBinLabel(11,"Vertex");
  fEvtCount->GetXaxis()->SetBinLabel(12,"Centrality");
  fEvtCount->GetXaxis()->SetBinLabel(13,"Pile up");
  fEvtCount->GetXaxis()->SetBinLabel(14,"Recenter");

  fQAList->Add(fEvtCount);

  ////////////////////////
  // Run Number Info
  ////////////////////////
  TString runNumList10h[90] = {
    "139510", "139507", "139505", "139503", "139465", "139438", "139437", "139360", "139329", "139328",
    "139314", "139310", "139309", "139173", "139107", "139105", "139038", "139037", "139036", "139029",
    "139028", "138872", "138871", "138870", "138837", "138732", "138730", "138666", "138662", "138653",
    "138652", "138638", "138624", "138621", "138583", "138582", "138578", "138534", "138469", "138442",
    "138439", "138438", "138396", "138364", "138275", "138225", "138201", "138197", "138192", "138190",
    "137848", "137844", "137752", "137751", "137724", "137722", "137718", "137704", "137693", "137692",
    "137691", "137686", "137685", "137639", "137638", "137608", "137595", "137549", "137546", "137544",
    "137541", "137539", "137531", "137530", "137443", "137441", "137440", "137439", "137434", "137432",
    "137431", "137430", "137243", "137236", "137235", "137232", "137231", "137230", "137162", "137161"};
  TString runNumList15o[138] = {
    "246994", "246991", "246989", "246984", "246982", "246948", "246945", "246928", "246871", "246870",
    "246867", "246865", "246864", "246859", "246858", "246851", "246847", "246846", "246845", "246844",
    "246810", "246809", "246808", "246807", "246805", "246804", "246766", "246765", "246763", "246760",
    "246759", "246758", "246757", "246751", "246750", "246434", "246431", "246424", "246392", "246391",
    "246276", "246275", "246272", "246271", "246225", "246222", "246217", "246185", "246182", "246181",
    "246180", "246178", "246153", "246152", "246151", "246148", "246115", "246113", "246089", "246087",
    "246053", "246052", "246049", "246048", "246042", "246037", "246036", "246012", "246003", "246001",
    "245963", "245954", "245952", "245949", "245923", "245833", "245831", "245829", "245793", "245785",
    "245775", "245766", "245759", "245752", "245731", "245729", "245705", "245702", "245692", "245683",
    "245554", "245545", "245544", "245543", "245542", "245540", "245535", "245507", "245505", "245504",
    "245501", "245497", "245496", "245454", "245453", "245450", "245446", "245441", "245411", "245410",
    "245409", "245407", "245401", "245397", "245396", "245353", "245349", "245347", "245346", "245345",
    "245343", "245259", "245233", "245232", "245231", "245152", "245151", "245146", "245145", "245068",
    "245066", "245064", "244983", "244982", "244980", "244975", "244918", "244917"};
  TString runNumList18q[125] = {
    "296623","296622","296621","296619","296618","296616","296615","296594","296553","296552",
    "296551","296550","296548","296547","296516","296512","296511","296510","296509","296472",
    "296433","296424","296423","296420","296419","296415","296414","296383","296381","296380",
    "296379","296378","296377","296376","296375","296312","296309","296304","296303","296280",
    "296279","296273","296270","296269","296247","296246","296244","296243","296242","296241",
    "296240","296198","296197","296196","296195","296194","296192","296191","296143","296142",
    "296135","296134","296133","296132","296123","296074","296066","296065","296063","296062",
    "296060","296016","295942","295941","295937","295936","295913","295910","295909","295861",
    "295860","295859","295856","295855","295854","295853","295831","295829","295826","295825",
    "295822","295819","295818","295816","295791","295788","295786","295763","295762","295759",
    "295758","295755","295754","295725","295723","295721","295719","295718","295717","295714",
    "295712","295676","295675","295673","295668","295667","295666","295615","295612","295611",
    "295610","295589","295588","295586","295585"};
  TString runNumList18r[89] = {
    "297595","297590","297588","297558","297544","297542","297541","297540","297537","297512",
    "297483","297479","297452","297451","297450","297446","297442","297441","297415","297414",
    "297413","297406","297405","297380","297379","297372","297367","297366","297363","297336",
    "297335","297333","297332","297317","297311","297310","297278","297222","297221","297218",
    "297196","297195","297193","297133","297132","297129","297128","297124","297123","297119",
    "297118","297117","297085","297035","297031","296966","296941","296938","296935","296934",
    "296932","296931","296930","296903","296900","296899","296894","296852","296851","296850",
    "296848","296839","296838","296836","296835","296799","296794","296793","296790","296787",
    "296786","296785","296784","296781","296752","296694","296693","296691","296690"};
  runNumList = new std::map<int,int>;
  if      (fPeriod.EqualTo("LHC10h")) for (int i = 0; i < 90; i++) runNumList->insert(std::pair<int,int>(runNumList10h[i].Atoi(),i+1));
  else if (fPeriod.EqualTo("LHC15o")) for (int i = 0; i <138; i++) runNumList->insert(std::pair<int,int>(runNumList15o[i].Atoi(),i+1));
  else if (fPeriod.EqualTo("LHC18q")) for (int i = 0; i <125; i++) runNumList->insert(std::pair<int,int>(runNumList18q[i].Atoi(),i+1));
  else if (fPeriod.EqualTo("LHC18r")) for (int i = 0; i < 89; i++) runNumList->insert(std::pair<int,int>(runNumList18r[i].Atoi(),i+1));
  else return;
  fHistRunNumBin = new TH1I("runNumBin","",(int)runNumList->size(),1,(int)runNumList->size()+1);
  std::map<int,int>::iterator iter;
  for (auto runNum : *runNumList) fHistRunNumBin->GetXaxis()->SetBinLabel(runNum.second, Form("%i",runNum.first));
  fQAList->Add(fHistRunNumBin);

  // Event-wise QA
  fHistCent[0] = new TH1D("fHistCentBfCut", "Dist. of Centrality Before Cut", 80, 0., 80.);
  fHistCent[1] = new TH1D("fHistCentAfCut", "Dist. of Centrality After Cut", 80, 0., 80.);
  fHistVz[0] = new TH1D("fHistVzBfCut", "Dist of Vz Before Cut", 200, -20., 20.);
  fHistVz[1] = new TH1D("fHistVzAfCut", "Dist of Vz After Cut", 200, -20., 20.);
  fQAList->Add(fHistCent[0]);
  fQAList->Add(fHistCent[1]);
  fQAList->Add(fHistVz[0]);
  fQAList->Add(fHistVz[1]);

  fProfileChanalMult[0] = new TProfile("fProfileChanalMult_RAW", ";Chanal;Energy", 64, 0., 64.);
  fProfileChanalMult[1] = new TProfile("fProfileChanalMult_BGE", ";Chanal;Energy", 64, 0., 64.);
  fQAList->Add(fProfileChanalMult[0]);
  fQAList->Add(fProfileChanalMult[1]);
  
  std::string name;
  for (int i = 0; i < 2; i++) {
    if (i == 0) name = "BfRC";
    if (i == 1) name = "AfRC";
    fProfileQ2xV0CCent[i] = new TProfile(Form("fProfileQ2xV0CCent_%s",name.c_str()),Form("fProfileQ2xV0CCent_%s;Cent;<Q2x>",name.c_str()),fCentBinNum,0.,80);
    fProfileQ2yV0CCent[i] = new TProfile(Form("fProfileQ2yV0CCent_%s",name.c_str()),Form("fProfileQ2yV0CCent_%s;Cent;<Q2y>",name.c_str()),fCentBinNum,0.,80);
    fProfileQ2xV0ACent[i] = new TProfile(Form("fProfileQ2xV0ACent_%s",name.c_str()),Form("fProfileQ2xV0ACent_%s;Cent;<Q2x>",name.c_str()),fCentBinNum,0.,80);
    fProfileQ2yV0ACent[i] = new TProfile(Form("fProfileQ2yV0ACent_%s",name.c_str()),Form("fProfileQ2yV0ACent_%s;Cent;<Q2y>",name.c_str()),fCentBinNum,0.,80);

    fProfileQ2xV0CVz[i] = new TProfile(Form("fProfileQ2xV0CVz_%s",name.c_str()),Form("fProfileQ2xV0CVz_%s;Vz;<Q2x>",name.c_str()),10,-10,10);
    fProfileQ2yV0CVz[i] = new TProfile(Form("fProfileQ2yV0CVz_%s",name.c_str()),Form("fProfileQ2yV0CVz_%s;Vz;<Q2y>",name.c_str()),10,-10,10);
    fProfileQ2xV0AVz[i] = new TProfile(Form("fProfileQ2xV0AVz_%s",name.c_str()),Form("fProfileQ2xV0AVz_%s;Vz;<Q2x>",name.c_str()),10,-10,10);
    fProfileQ2yV0AVz[i] = new TProfile(Form("fProfileQ2yV0AVz_%s",name.c_str()),Form("fProfileQ2yV0AVz_%s;Vz;<Q2y>",name.c_str()),10,-10,10);

    fQAList->Add(fProfileQ2xV0CCent[i]);
    fQAList->Add(fProfileQ2yV0CCent[i]);
    fQAList->Add(fProfileQ2xV0ACent[i]);
    fQAList->Add(fProfileQ2yV0ACent[i]);

    fQAList->Add(fProfileQ2xV0CVz[i]);
    fQAList->Add(fProfileQ2yV0CVz[i]);
    fQAList->Add(fProfileQ2xV0AVz[i]);
    fQAList->Add(fProfileQ2yV0AVz[i]);
  }
  

  fResultsList = new TList();
  fResultsList -> SetName("fResultsList");
  fResultsList -> SetOwner(kTRUE);

  //////////////////////////Recenter file/////////////////////////////////
  for (std::map<int,int>::reverse_iterator r_it = runNumList->rbegin(); r_it != runNumList->rend(); ++r_it) {
    p2D_V0C_meanQ2x_cent_vz[r_it->second - 1] = new TProfile2D(Form("p2D_V0C_meanQ2x_cent_vz_run%i",r_it->first),Form("p2D_V0C_meanQ2x_cent_vz_run%i",r_it->first),fCentBinNum,0,80,fVzBinNum,-10,10);
    p2D_V0C_meanQ2y_cent_vz[r_it->second - 1] = new TProfile2D(Form("p2D_V0C_meanQ2y_cent_vz_run%i",r_it->first),Form("p2D_V0C_meanQ2y_cent_vz_run%i",r_it->first),fCentBinNum,0,80,fVzBinNum,-10,10);
    p2D_V0A_meanQ2x_cent_vz[r_it->second - 1] = new TProfile2D(Form("p2D_V0A_meanQ2x_cent_vz_run%i",r_it->first),Form("p2D_V0A_meanQ2x_cent_vz_run%i",r_it->first),fCentBinNum,0,80,fVzBinNum,-10,10);
    p2D_V0A_meanQ2y_cent_vz[r_it->second - 1] = new TProfile2D(Form("p2D_V0A_meanQ2y_cent_vz_run%i",r_it->first),Form("p2D_V0A_meanQ2y_cent_vz_run%i",r_it->first),fCentBinNum,0,80,fVzBinNum,-10,10);
    fResultsList->Add(p2D_V0C_meanQ2x_cent_vz[r_it->second - 1]);
    fResultsList->Add(p2D_V0C_meanQ2y_cent_vz[r_it->second - 1]);
    fResultsList->Add(p2D_V0A_meanQ2x_cent_vz[r_it->second - 1]);
    fResultsList->Add(p2D_V0A_meanQ2y_cent_vz[r_it->second - 1]);
  }
  
  PostData(1,fQAList);
  PostData(2,fResultsList);
  if (fDebug) Printf("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskCVEVZEROCalib::UserExec(Option_t *)
{
  if (fDebug) Printf("===============================We are in UserExec!================================");
  fEvtCount->Fill(1);
  //----------------------------
  // Handle
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else fEvtCount->Fill(2);

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else fEvtCount->Fill(3);

  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else fEvtCount->Fill(4);

  fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else fEvtCount->Fill(6);

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel) {
      AliError(Form("%s: Could not get AliMultSelection", GetName()));
    } else fEvtCount->Fill(7);
    if (!manager || !handler || !fAOD || !fUtils || !fMultSel) return;
  }

  if (!manager || !handler || !fAOD || !fUtils) return;
  if (fDebug) Printf("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  UInt_t mask = handler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB"))
  isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7"))
  isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  if (isTrigselected == false) return;
  fEvtCount->Fill(8);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Run Number
  //----------------------------
  fRunNum = fAOD->GetRunNumber();
  fEvtCount->Fill(9);
  if (fRunNum != fOldRunNum) {
     // Load the run dependent calibration hist
      if (!LoadCalibHistForThisRun()) return;
      fRunNumBin = runNumList->at(fRunNum);
      fOldRunNum = fRunNum;
      if (fRunNumBin < 0) return;
  }
  fHistRunNumBin->Fill(fRunNumBin);
  fEvtCount->Fill(10);
  if (fDebug) Printf("read in done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  fVtx -> GetXYZ(fVertex);
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (fabs(fVertex[0])<1e-6 || fabs(fVertex[1])<1e-6 || fabs(fVertex[2])<1e-6) return;
  double dz = fVertex[2] - fAOD->GetPrimaryVertexSPD()->GetZ();
  fHistVz[0]->Fill(fVertex[2]);
  if (fabs(fVertex[2]) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  // fEventCuts.SetCentralityEstimators("V0M","CL1");
  // if (!fEventCuts->AcceptEvent(fAOD) ) return;
  if (fPeriod.EqualTo("LHC10h")) if (fabs(dz)>0.5) return;
  if (fPeriod.EqualTo("LHC15o")) {
      double covTrc[6],covSPD[6];
      fVtx->GetCovarianceMatrix(covTrc);
      fAOD->GetPrimaryVertexSPD()->GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (fabs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return;
  }
  fHistVz[1]->Fill(fVertex[2]);
  for (int i = 0; i < 20; ++i) {
      if (fVertex[2] > -10+i*1 && fVertex[2] < -10+(i+1)*1) {fVzBin = i; break;}
  }
  if (fVzBin<0) return;
  fEvtCount->Fill(11);
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  if (fPeriod.EqualTo("LHC10h")) {
    fCentV0M  = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    fCentTRK  = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    fCentSPD0 = fAOD->GetCentrality()->GetCentralityPercentile("CL0");
    fCentSPD1 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
  } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fCentV0M  = fMultSel->GetMultiplicityPercentile("V0M");
    fCentTRK  = fMultSel->GetMultiplicityPercentile("TRK");
    fCentSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
    fCentSPD1 = fMultSel->GetMultiplicityPercentile("CL1");
  } else return;
  //we use centV0M as the default centrality
  fCent = fCentV0M;
  if (fCent < 0 || fCent >= 80) return;
  // cent bin
  fCentBin = (int)fCent/10;
  fHistCent[0]->Fill(fCent);
  fEvtCount->Fill(12);
  if (fDebug) Printf("centrality done!");

  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC10h")) if (!RemovalForRun1()) return;
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    // hMultCentQA[0]->Fill(fCent, fAOD->GetNumberOfTracks()); // raw Trk Multi Vs Cent(V0M)
    // if (PileUpMultiVertex(fAOD)) return;
    // if (!RejectEvtMultComp(fAOD)) return;
    // hMultCentQA[1]->Fill(fCent, fAOD->GetNumberOfTracks()); // Mult_Cent QA
    // if (!AODPileupCheck (fAOD)) return;
    if (!RejectEvtTFFit()) return; // 15o_pass2
  }
  fHistCent[1]->Fill(fCent);
  fEvtCount->Fill(13);
  if (fDebug) Printf("pile-up done!");
  //----------------------------
  // VZERO Plane
  //----------------------------

  if(!RecenterVZEROQVector()) return;
  fEvtCount->Fill(14);
  if (fDebug) Printf("recenter done!");

  //------------------
  // Post output data.
  //------------------
  PostData(1,fQAList);
  PostData(2,fResultsList);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::RecenterVZEROQVector()
{
  double multMean[64] = {0};
  double Q2xV0C = 0.;
  double Q2yV0C = 0.;
  double multV0C = 0.;

  double Q2xV0A = 0.;
  double Q2yV0A = 0.;
  double multV0A = 0.;

  //Load the GE histograms
  if (fPeriod.EqualTo("LHC10h") ) 
  for (int iCh = 0; iCh < 64; ++iCh) multMean[iCh] = hMultV0Read->GetBinContent(iCh+1, fRunNumBin+1);
  if (fPeriod.EqualTo("LHC15o"))  
  for (int iCh = 0; iCh < 64; ++iCh) multMean[iCh] = hMultV0->GetBinContent(iCh+1);

  //Loop Over VZERO Channels
  //Gain Equalization
  for (int iCh = 0; iCh < 64; ++iCh) {
    double phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    double multCh = 0.;

    if      (fPeriod.EqualTo("LHC10h")) multCh = fAOD->GetVZEROEqMultiplicity(iCh);
    else if (fPeriod.EqualTo("LHC15o")|| 
             fPeriod.EqualTo("LHC18q")||
             fPeriod.EqualTo("LHC18r")) multCh = fAOD->GetVZEROData()->GetMultiplicity(iCh);
    else return false;
  
    if (iCh < 32) { // C-side
      if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC15o")) {
        fProfileChanalMult[0] -> Fill(iCh+0.5, multCh);
        if      (iCh <  8)              multCh = multCh/multMean[iCh] * multMean[0];
        else if (iCh >= 8  && iCh < 16) multCh = multCh/multMean[iCh] * multMean[8];
        else if (iCh >= 16 && iCh < 24) multCh = multCh/multMean[iCh] * multMean[16];
        else if (iCh >= 24 && iCh < 32) multCh = multCh/multMean[iCh] * multMean[24];
        fProfileChanalMult[1] -> Fill(iCh+0.5, multCh);
      }
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
        fProfileChanalMult[0] -> Fill(iCh+0.5, multCh);
        int    ibinV0     = fHCorrectV0ChWeghts->FindBin(fVertex[2],iCh);
        double gainFactor = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
        multCh *= gainFactor;
        fProfileChanalMult[1] -> Fill(iCh+0.5, multCh);
      }

      if (multCh < 1.e-6) return false;

      Q2xV0C  += multCh * TMath::Cos(2 * phi);
      Q2yV0C  += multCh * TMath::Sin(2 * phi);
      multV0C += multCh;

    } else if (iCh >= 32 && iCh < 64) { // A-side
      fProfileChanalMult[0] -> Fill(iCh+0.5, multCh);
      if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC15o")) {
        if      (iCh >= 32 && iCh < 40) multCh = multCh/multMean[iCh] * multMean[0];
        else if (iCh >= 40 && iCh < 48) multCh = multCh/multMean[iCh] * multMean[8];
        else if (iCh >= 48 && iCh < 56) multCh = multCh/multMean[iCh] * multMean[16];
        else if (iCh >= 56 && iCh < 64) multCh = multCh/multMean[iCh] * multMean[24];
        fProfileChanalMult[1] -> Fill(iCh+0.5, multCh);
      }
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
        fProfileChanalMult[0] -> Fill(iCh+0.5, multCh);
        int    ibinV0     = fHCorrectV0ChWeghts->FindBin(fVertex[2],iCh);
        double gainFactor = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
        multCh *= gainFactor;
        fProfileChanalMult[1] -> Fill(iCh+0.5, multCh);
      }

      if (multCh < 1.e-6) return false;

      Q2xV0A  += multCh * TMath::Cos(2 * phi);
      Q2yV0A  += multCh * TMath::Sin(2 * phi);
      multV0A += multCh;
    }
  }

  Q2xV0C /= multV0C;
  Q2yV0C /= multV0C;

  Q2xV0A /= multV0A;
  Q2yV0A /= multV0A;

  if(runNumList->find(fRunNum) == runNumList->end()) return false;
  int rumNumbin = runNumList->at(fRunNum) - 1;
  p2D_V0C_meanQ2x_cent_vz[rumNumbin] -> Fill(fCentSPD0, fVertex[2], Q2xV0C);
  p2D_V0C_meanQ2y_cent_vz[rumNumbin] -> Fill(fCentSPD0, fVertex[2], Q2yV0C);
  p2D_V0A_meanQ2x_cent_vz[rumNumbin] -> Fill(fCentSPD0, fVertex[2], Q2xV0A);
  p2D_V0A_meanQ2y_cent_vz[rumNumbin] -> Fill(fCentSPD0, fVertex[2], Q2yV0A);

  fProfileQ2xV0CCent[0]->Fill(fCentSPD0,Q2xV0C);
  fProfileQ2yV0CCent[0]->Fill(fCentSPD0,Q2yV0C);
  fProfileQ2xV0ACent[0]->Fill(fCentSPD0,Q2xV0A);
  fProfileQ2yV0ACent[0]->Fill(fCentSPD0,Q2yV0A);
  fProfileQ2xV0CVz[0]->Fill(fVertex[2],Q2xV0C);
  fProfileQ2yV0CVz[0]->Fill(fVertex[2],Q2yV0C);
  fProfileQ2xV0AVz[0]->Fill(fVertex[2],Q2xV0A);
  fProfileQ2yV0AVz[0]->Fill(fVertex[2],Q2yV0A);
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::LoadCalibHistForThisRun()
{
  if (fPeriod.EqualTo("LHC15o")) {
    // 15o VZERO Calibration Histograms
    hMultV0 -> Reset();
    for (int i = 0; i < 2; i++) {
      hQx2mV0[i] -> Reset();
      hQy2mV0[i] -> Reset();
    }
    hMultV0    = ((TH1D*) contMult ->GetObject(fRunNum));
    hQx2mV0[0] = ((TH1D*) contQxncm->GetObject(fRunNum));
    hQy2mV0[0] = ((TH1D*) contQyncm->GetObject(fRunNum));
    hQx2mV0[1] = ((TH1D*) contQxnam->GetObject(fRunNum));
    hQy2mV0[1] = ((TH1D*) contQynam->GetObject(fRunNum));
    if (!hMultV0)    return false;
    if (!hQx2mV0[0]) return false;
    if (!hQy2mV0[0]) return false;
    if (!hQx2mV0[1]) return false;
    if (!hQy2mV0[1]) return false;
  }

  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    //18q/r VZERO
    for (int i = 0; i < 2; i++) {
      hQx2mV0[i] ->Reset();
      hQy2mV0[i] ->Reset();
    }
    hQx2mV0[0] = ((TH1D*) contQxncm->GetObject(fRunNum));
    hQy2mV0[0] = ((TH1D*) contQyncm->GetObject(fRunNum));
    hQx2mV0[1] = ((TH1D*) contQxnam->GetObject(fRunNum));
    hQy2mV0[1] = ((TH1D*) contQynam->GetObject(fRunNum));
    for (int i = 0; i < 2; i++) {
      if (!hQx2mV0[i]) return false;
      if (!hQy2mV0[i]) return false;
    }
    fHCorrectV0ChWeghts -> Reset();
    fHCorrectV0ChWeghts = (TH2F *) fListVZEROCalib->FindObject(Form("hWgtV0ChannelsvsVzRun%d",fRunNum));
    if (!fHCorrectV0ChWeghts) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::RemovalForRun1()
{
  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  fUtils->SetUseMVPlpSelection(true);
  // fUtils->SetMinPlpContribMV(5);
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
  // bool isPileup = fUtils->IsPileUpMV(fAOD); // pp, p-Pb
  if (isPileup) return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::RejectEvtMultComp() // 15o_pass1, old pile-up
{
   // TPC cluster cut
    Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks(); // multESD
    const Int_t nTrk = fAOD->GetNumberOfTracks();
    Int_t multTPC=0; // FB128 + Common Track Cuts
    Int_t multTPCFE=0; // FB1 + Common Track Cuts + chi2 > 0.2
    Int_t multGlobal=0; // FB16 + Common Track Cuts + Dca Cuts
    for (Int_t it1 = 0; it1 < nTrk; it1++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it1);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(128)) multTPC++;
      double Eta  = aodTrk->Eta();
      double Pt    = aodTrk->Pt();
      // double Phi  = aodTrk->Phi();
      if (Pt<0.2 || Pt>5.0 || TMath::Abs(Eta)>0.8 || aodTrk->GetTPCNcls()<70. || aodTrk->GetTPCsignal()<10.0) continue;
      if (aodTrk->TestFilterBit(1) && aodTrk->Chi2perNDF()>0.2)  multTPCFE++;
      if (!aodTrk->TestFilterBit(16) || aodTrk->Chi2perNDF()<0.1)   continue;
      Double_t dca[2] = {-99., -99.};
      Double_t cov[3] = {-99., -99., -99.};
      Double_t magField = fAOD->GetMagneticField();
      if (magField!=0) {
        if (aodTrk->PropagateToDCA(fAOD->GetPrimaryVertex(), magField, 100., dca, cov) && TMath::Abs(dca[0]) < 0.3 && TMath::Abs(dca[1]) < 0.3) multGlobal++;
      }
    }
    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::RejectEvtTFFit()
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

  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::RejectEvtTPCITSfb32TOF ()
{
    //TOD+FB32 pile-up removal
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
    Int_t multTrk=0;
    Int_t multTrkTOF=0;
    int nTrk = fAOD->GetNumberOfTracks();
    for (Int_t it2 = 0; it2 < nTrk; it2++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it2);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(32)) {
        multTrk++;
        if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10 && aodTrk->GetTOFsignal() >= 12000 && aodTrk->GetTOFsignal() <= 25000) multTrkTOF++;
        else return false;
      }
    }
    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::AODPileupCheck()
{
  Int_t isPileup = fAOD->IsPileupFromSPD(3);
  if (isPileup !=0 && fPeriod.EqualTo("LHC16t")) return false; // LHC16t : pPb
  if (fAOD->IsIncompleteDAQ()) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fPeriod.EqualTo("LHC15o")) {
    if (!fMultSel->GetThisEventIsNotPileup()) return false;
    if (!fMultSel->GetThisEventIsNotPileupMV()) return false;
    if (!fMultSel->GetThisEventIsNotPileupInMultBins()) return false;
    if (!fMultSel->GetThisEventHasNoInconsistentVertices()) return false;
    if (!fMultSel->GetThisEventPassesTrackletVsCluster()) return false;
    if (!fMultSel->GetThisEventIsNotIncompleteDAQ()) return false;
    if (!fMultSel->GetThisEventHasGoodVertex2016()) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEVZEROCalib::PileUpMultiVertex()
{
  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if (!(nPlp=fAOD->GetNumberOfPileupVerticesTracks()))
  return false;

  vtPrm = fAOD->GetPrimaryVertex();
  if (vtPrm == fAOD->GetPrimaryVertexSPD())
  return true;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)fAOD->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;

    return true; // pile-up: well separated vertices
  }
  return false;
}

//---------------------------------------------------

double AliAnalysisTaskCVEVZEROCalib::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
    // calculate sqrt of weighted distance to other vertex
    if (!v0 || !v1) {
        printf("One of vertices is not valid\n");
        return 0;
    }
    static TMatrixDSym vVb(3);
    double dist = -1;
    double dx = v0->GetX()-v1->GetX();
    double dy = v0->GetY()-v1->GetY();
    double dz = v0->GetZ()-v1->GetZ();
    double cov0[6],cov1[6];
    v0->GetCovarianceMatrix(cov0);
    v1->GetCovarianceMatrix(cov1);
    vVb(0,0) = cov0[0]+cov1[0];
    vVb(1,1) = cov0[2]+cov1[2];
    vVb(2,2) = cov0[5]+cov1[5];
    vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
    vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
    vVb.InvertFast();
    if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
    dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
        +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
    return dist>0 ? TMath::Sqrt(dist) : -1;
}

