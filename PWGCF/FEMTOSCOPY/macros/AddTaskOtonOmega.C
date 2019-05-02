#include "TROOT.h"
#include "TSystem.h"
AliAnalysisTaskSE* AddTaskOtonOmega(	bool isMC = false,     
                                     	bool isESD = false,
                                     	TString CentEst = "kInt7",
                                     	bool CascadeTreeFlag = false,
                                     	bool OmegaTreeFlag = false,
                                     	Int_t CutFlag = 0,
                                     	bool GetConfigFromAlien = true,
                                    	TString cFileName = "ConfigOtonOmega.C", 
					bool DCAPlots = false,
					bool DeltaEtaDeltaPhiCut = false,
					bool RunNumberQA = false
)
{


//Start with fixed definitions from AddTaskFemtoDream:
bool notpp = true;  //1
bool fineBinning = true;  //2
bool CPAPlots = false;  //4
bool MomReso = false;  //5
bool etaPhiPlotsAtTPCRadii = false;  //6
bool CombSigma = false;  //7
bool PileUpRej = true;  //8
bool mTkTPlot = true;  //9
bool kTCentPlot = false;  //10
bool MultvsCentPlot = false;  //11
bool dPhidEtaPlots = false;  //12
bool eventMixing = true;  //13
bool phiSpin = true;  //14
bool stravinskyPhiSpin = true;  //15
bool ContributionSplitting = false;  //16
bool ContributionSplittingDaug = false;  //17
int FilterBit = 128;  //19
bool InvMassPairs = false;  //20
int SphericityRange = 0;  //22



    //Get Config File:
    TString configBasePath= "$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/";
    if(GetConfigFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/o/ovazquez/%s .",cFileName.Data()))) ){
        configBasePath=Form("%s/",gSystem->pwd());
    } else {
        cFileName = "ConfigOtonOmega.C"; //if not getting it from alien, take the standard one
    }
    TString configFilePath(configBasePath+cFileName);
    std::cout << "Configpath: " << configFilePath << std::endl;
    TString configFunction(cFileName(0,cFileName.Length() - 2));
    if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFunction.Data())) gROOT->LoadMacro(configFilePath.Data());





  // 1    2     3     4     5     6     7    8    9      10   11     12   13    14    15    16   17
  //true,true,false,false,false,false,false,true,false,false,true,false,true,false,false,false,true
  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  if (!(AliPIDResponse*) mgr->GetTask("PIDResponseTask")) {
    if (isMC) {
      // IMPORTANT - SET WHEN USING DIFFERENT PASS
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
          gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/"
                                     "AddTaskPIDResponse.C (kTRUE, kTRUE, "
                                     "kTRUE, \"1\")"));
    } else {
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
          gInterpreter->ExecuteMacro(
              "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C)"));
    }
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetMultVsCentPlots(MultvsCentPlot);
  if (SphericityRange != 0) {
    if (SphericityRange == 1) {
      evtCuts->SetSphericityCuts(0., 0.3);
    } else if (SphericityRange == 2) {
      evtCuts->SetSphericityCuts(0.3, 0.7);
    } else if (SphericityRange == 3) {
      evtCuts->SetSphericityCuts(0.7, 1.0);
    } else if (SphericityRange == 4) {
      evtCuts->SetSphericityCuts(0.7, 0.8);
    } else if (SphericityRange == 5) {
      evtCuts->SetSphericityCuts(0.8, 0.9);
    } else if (SphericityRange == 6) {
      evtCuts->SetSphericityCuts(0.9, 1.0);
    } else {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      std::cout << "Warning: Unkown Sphericity Type\n";
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

    }
  }


  //Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, DCAPlots, CombSigma, ContributionSplitting);
  if (isESD) {
    TrackCuts->SetCheckFilterBit(false);
    TrackCuts->SetCheckESDFiltering(true);
    TrackCuts->SetDCAReCalculation(false);
  } else {
    TrackCuts->SetFilterBit(FilterBit);
  }
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, DCAPlots, CombSigma,
                                             ContributionSplitting);
  if (isESD) {
    AntiTrackCuts->SetCheckFilterBit(false);
    AntiTrackCuts->SetCheckESDFiltering(true);
    AntiTrackCuts->SetDCAReCalculation(false);
  } else {
    AntiTrackCuts->SetFilterBit(FilterBit);
  }
  AntiTrackCuts->SetCutCharge(-1);


  //obsolete v0 cuts:
  /*
  AliFemtoDreamv0Cuts *v0Cuts;
  AliFemtoDreamv0Cuts *Antiv0Cuts;
  v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, CPAPlots,
                                           ContributionSplitting);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, PileUpRej, false);
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, PileUpRej, false);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda
  Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, CPAPlots,
                                               ContributionSplitting);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, PileUpRej, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, PileUpRej, false);
  NegAntiv0Daug->SetCutCharge(-1);
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda
  */


  AliOtonOmegaCascadeCuts *CascadeCuts;
  CascadeCuts = AliOtonOmegaCascadeCuts::XiCuts(isMC, ContributionSplitting, CutFlag);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, PileUpRej, false);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, PileUpRej, false);
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(isMC, PileUpRej, false);
  //AliOtonOmegaTrackCuts *XiBachCuts = AliOtonOmegaTrackCuts::OmegaBachKaonCuts(isMC, PileUpRej, false);
  //CascadeCuts->SetPDGCodeCasc(3312);
  CascadeCuts->SetPDGCodeCasc(3334);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  //CascadeCuts->SetPDGCodeBach(-211);
  CascadeCuts->SetPDGCodeBach(-321);


  AliOtonOmegaCascadeCuts *AntiCascadeCuts;
  AntiCascadeCuts = AliOtonOmegaCascadeCuts::XiCuts( isMC, ContributionSplitting, CutFlag);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, PileUpRej, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, PileUpRej, false);
  AntiXiPosCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiXiBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(isMC, PileUpRej, false);
  //AliOtonOmegaTrackCuts *AntiXiBachCuts = AliOtonOmegaTrackCuts::OmegaBachKaonCuts(isMC, PileUpRej, false);
  AntiXiBachCuts->SetCutCharge(1);
  //AntiCascadeCuts->SetPDGCodeCasc(-3312);
  AntiCascadeCuts->SetPDGCodeCasc(-3334);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  //AntiCascadeCuts->SetPDGCodeBach(-211);
  AntiCascadeCuts->SetPDGCodeBach(-321);


  AliOtonOmegaCascadeCuts *CascadeOmegaCuts;
  CascadeOmegaCuts = AliOtonOmegaCascadeCuts::OmegaCuts(isMC,ContributionSplitting, CutFlag);
  CascadeOmegaCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *OmegaNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC,PileUpRej,false);
  AliFemtoDreamTrackCuts *OmegaPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC,PileUpRej,false);
  AliFemtoDreamTrackCuts *OmegaBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(isMC,PileUpRej,false);
  //AliOtonOmegaTrackCuts *OmegaBachCuts = AliOtonOmegaTrackCuts::OmegaBachKaonCuts(isMC,PileUpRej,false);
  CascadeOmegaCuts->SetPDGCodeCasc(3334);
  CascadeOmegaCuts->SetPDGCodev0(3122);
  CascadeOmegaCuts->SetPDGCodePosDaug(2212);
  CascadeOmegaCuts->SetPDGCodeNegDaug(-211);
  CascadeOmegaCuts->SetPDGCodeBach(-321);


  AliOtonOmegaCascadeCuts *AntiCascadeOmegaCuts;
  AntiCascadeOmegaCuts = AliOtonOmegaCascadeCuts::OmegaCuts(isMC,ContributionSplitting, CutFlag);
  AntiCascadeOmegaCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiOmegaNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC,PileUpRej,false);
  AntiOmegaNegCuts->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *AntiOmegaPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC,PileUpRej,false);
  AntiOmegaPosCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiOmegaBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(isMC,PileUpRej,false);
  //AliOtonOmegaTrackCuts *AntiOmegaBachCuts = AliOtonOmegaTrackCuts::OmegaBachKaonCuts(isMC,PileUpRej,false);
  AntiOmegaBachCuts->SetCutCharge(1);
  AntiCascadeOmegaCuts->SetPDGCodeCasc(-3334);
  AntiCascadeOmegaCuts->SetPDGCodev0(-3122);
  AntiCascadeOmegaCuts->SetPDGCodePosDaug(211);
  AntiCascadeOmegaCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeOmegaCuts->SetPDGCodeBach(-321);


  //Pass cuts to config file (to be further setup of modified there)_
  ConfigOtonOmega(CascadeCuts,XiPosCuts,XiNegCuts,XiBachCuts,AntiCascadeCuts,AntiXiPosCuts,AntiXiNegCuts,AntiXiBachCuts,CascadeOmegaCuts,OmegaPosCuts,OmegaNegCuts,OmegaBachCuts,AntiCascadeOmegaCuts,AntiOmegaPosCuts,AntiOmegaNegCuts,AntiOmegaBachCuts,TrackCuts,AntiTrackCuts);


  //finalize cuts set:
  if(RunNumberQA) {
   CascadeCuts->SetRunNumberQA(252234, 294926);
   AntiCascadeCuts->SetRunNumberQA(252234, 294926);
   CascadeOmegaCuts->SetRunNumberQA(252234, 294926);
   AntiCascadeOmegaCuts->SetRunNumberQA(252234, 294926);
  }
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  CascadeOmegaCuts->Setv0Negcuts(OmegaNegCuts);
  CascadeOmegaCuts->Setv0PosCuts(OmegaPosCuts);
  CascadeOmegaCuts->SetBachCuts(OmegaBachCuts);
  AntiCascadeOmegaCuts->Setv0Negcuts(AntiOmegaNegCuts);
  AntiCascadeOmegaCuts->Setv0PosCuts(AntiOmegaPosCuts);
  AntiCascadeOmegaCuts->SetBachCuts(AntiOmegaBachCuts);

  //skip v0s
  //if (RunNumberQA) {
    //    v0Cuts->SetRunNumberQA(265309, 267167);
    //    Antiv0Cuts->SetRunNumberQA(265309, 267167);
  //}

  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
//reorder to get Xi and Omega, this order is preserved afterwards in the code,
//for example in AliOtonOmegaAnalysis when doing the paircleaner
//  PDGParticles.push_back(3122);
//  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3334);
  PDGParticles.push_back(3334);
  //std::vector<double> ZVtxBins = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(-8);
  ZVtxBins.push_back(-6);
  ZVtxBins.push_back(-4);
  ZVtxBins.push_back(-2);
  ZVtxBins.push_back(0);
  ZVtxBins.push_back(2);
  ZVtxBins.push_back(4);
  ZVtxBins.push_back(6);
  ZVtxBins.push_back(8);
  ZVtxBins.push_back(10);
  std::vector<int> NBins;
  if (fineBinning) {
    NBins.push_back(750);  // p p
    NBins.push_back(750);  // p barp
    NBins.push_back(750);  // p Lambda
    NBins.push_back(750);  // p barLambda
    NBins.push_back(750);  // p Xi
    NBins.push_back(750);  // p barXi
    NBins.push_back(750);  // barp barp
    NBins.push_back(750);  // barp Lambda
    NBins.push_back(750);  // barp barLambda
    NBins.push_back(750);  // barp Xi
    NBins.push_back(750);  // barp barXi
    NBins.push_back(750);  // Lambda Lambda
    NBins.push_back(750);  // Lambda barLambda
    NBins.push_back(750);  // Lambda Xi
    NBins.push_back(750);  // Lambda barXi
    NBins.push_back(750);  // barLambda barLambda
    NBins.push_back(750);  // barLambda Xi
    NBins.push_back(750);  // barLambda barXi
    NBins.push_back(750);  // Xi Xi
    NBins.push_back(750);  // Xi barXi
    NBins.push_back(750);  // barXi barXi
  } else {  //standard binning Run1
    NBins.push_back(750);  // p p
    NBins.push_back(750);  // p barp
    NBins.push_back(150);  // p Lambda
    NBins.push_back(150);  // p barLambda
    NBins.push_back(150);  // p Xi
    NBins.push_back(150);  // p barXi
    NBins.push_back(750);  // barp barp
    NBins.push_back(150);  // barp Lambda
    NBins.push_back(150);  // barp barLambda
    NBins.push_back(150);  // barp Xi
    NBins.push_back(150);  // barp barXi
    NBins.push_back(150);  // Lambda Lambda
    NBins.push_back(150);  // Lambda barLambda
    NBins.push_back(150);  // Lambda Xi
    NBins.push_back(150);  // Lambda barXi
    NBins.push_back(150);  // barLambda barLambda
    NBins.push_back(150);  // barLambda Xi
    NBins.push_back(150);  // barLambda barXi
    NBins.push_back(150);  // Xi Xi
    NBins.push_back(150);  // Xi barXi
    NBins.push_back(150);  // barXi barXi
  }
  std::vector<float> kMin;
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");
  if (notpp) {
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(44);
    MultBins.push_back(48);
    MultBins.push_back(52);
    MultBins.push_back(56);
    MultBins.push_back(60);
    MultBins.push_back(64);
    MultBins.push_back(68);
    MultBins.push_back(72);
    MultBins.push_back(76);
    MultBins.push_back(80);
    MultBins.push_back(84);
    MultBins.push_back(88);
    MultBins.push_back(92);
    MultBins.push_back(96);
    MultBins.push_back(100);
    config->SetMultBins(MultBins);
  } else {
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(60);
    MultBins.push_back(80);
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);
  if (notpp)
    config->SetCentBinning(true);
  config->SetkTBinning(mTkTPlot);
  config->SetmTBinning(mTkTPlot);
  config->SetkTCentralityBinning(kTCentPlot);
  config->SetInvMassPairs(InvMassPairs);
  if (kTCentPlot) {
    std::vector<float> centBins;
    centBins.push_back(20);
    centBins.push_back(40);
    centBins.push_back(90);
    config->SetCentBins(centBins);
  }
  config->SetZBins(ZVtxBins);
  if (MomReso) {
    if (isMC) {
      config->SetMomentumResolution(true);
    } else {
      std::cout
          << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
    }
  }
  if (etaPhiPlotsAtTPCRadii) {
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout
          << "You are trying to request the Eta Phi Plots without MC Info; fix it wont work! \n";
    }
  }
  //New delta eta deta phi:
  if (DeltaEtaDeltaPhiCut) {
   std::vector<bool> CPR;
   CPR.push_back(true);  // p p
   CPR.push_back(true);  // p barp
   CPR.push_back(true);  // p Lambda
   CPR.push_back(true);  // p barLambda
   CPR.push_back(true);  // p Xi
   CPR.push_back(true);  // p barXi
   CPR.push_back(true);  // barp barp
   CPR.push_back(true);  // barp Lambda
   CPR.push_back(true);  // barp barLambda
   CPR.push_back(true);  // barp Xi
   CPR.push_back(true);  // barp barXi
   CPR.push_back(true);  // Lambda Lambda
   CPR.push_back(true);  // Lambda barLambda
   CPR.push_back(true);  // Lambda Xi
   CPR.push_back(true);  // Lambda barXi
   CPR.push_back(true);  // barLambda barLambda
   CPR.push_back(true);  // barLambda Xi
   CPR.push_back(true);  // barLambda barXi
   CPR.push_back(true);  // Xi Xi
   CPR.push_back(true);  // Xi barXi
   CPR.push_back(true);  // barXi barXi
   config->SetClosePairRejection(CPR);
   config->SetDeltaEtaMax(0.01);
   config->SetDeltaPhiMax(0.01);
  }
  config->SetdPhidEtaPlots(dPhidEtaPlots);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetSpinningDepth(10);
  config->SetUseEventMixing(eventMixing);
  config->SetUsePhiSpinning(phiSpin);
  config->SetUseStravinskyMethod(stravinskyPhiSpin);
  config->SetMinimalBookingME(false);
  config->SetMinimalBookingSample(true);
  if (!notpp) {
    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kSPD);
  } else {
    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  }
  AliAnalysisTaskOtonOmega *task = new AliAnalysisTaskOtonOmega(
      "FemtoDreamDefault", isESD, isMC, CascadeTreeFlag, OmegaTreeFlag);
  if (CentEst == "kInt7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
    task->SetMVPileUp(kTRUE);
  } else if (CentEst == "kMB") {
    task->SelectCollisionCandidates(AliVEvent::kMB);
    std::cout << "Added kMB Trigger \n";
    task->SetMVPileUp(kFALSE);
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
    task->SetMVPileUp(kFALSE);
  } else {
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "Centrality Estimator not set, fix it else your Results will be empty!"
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
  }
  //	task->SetDebugLevel(0);
  task->SetEvtCutQA(true);
  task->SetTrackBufferSize(2000);
  task->SetEventCuts(evtCuts);
  task->SetTrackCuts(TrackCuts);
  task->SetAntiTrackCuts(AntiTrackCuts);
//again skip v0s
//  task->Setv0Cuts(v0Cuts);
//  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetCascadeCuts(CascadeCuts);
  task->SetAntiCascadeCuts(AntiCascadeCuts);
  task->SetCascadeOmegaCuts(CascadeOmegaCuts);
  task->SetAntiCascadeOmegaCuts(AntiCascadeOmegaCuts);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString addon = "";    
  if (CentEst == "kInt7") {
    addon += "MB";
  } else if (CentEst == "kHM") {
    addon += "HM";
  }
  if (SphericityRange != 0 ) {
    addon += "_Sphericity_";
    addon += SphericityRange;
    addon += "_";
  }

  if (cFileName == "ConfigOtonOmega.C"){
    std::cout << "No CONTAINERaddon from configfile " << std::endl;
  }else{
    addon += "_";
    for(int ii=15;ii<=22;ii++){ addon += cFileName(ii);};
    std::cout << "CONTAINERaddon from configfile " << cFileName(15,22).Data() << std::endl;
  }

  std::cout << "CONTAINTER NAME: " << addon.Data() << std::endl;
  TString QAName = Form("%sQA", addon.Data());
  coutputQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      QAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  AliAnalysisDataContainer *coutputEvtCuts;
  TString EvtCutsName = Form("%sEvtCuts", addon.Data());
  coutputEvtCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      EvtCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *couputTrkCuts;
  TString TrackCutsName = Form("%sTrackCuts", addon.Data());
  couputTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      TrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  AliAnalysisDataContainer *coutputAntiTrkCuts;
  TString AntiTrackCutsName = Form("%sAntiTrackCuts", addon.Data());
  coutputAntiTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiTrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

//skip v0s
/*
  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts", addon.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts", addon.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);
*/

  AliAnalysisDataContainer *coutputCascadeCuts;
  TString CascadeCutsName = Form("%sCascadeCuts", addon.Data());
  coutputCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputCascadeCuts);

  AliAnalysisDataContainer *coutputAntiCascadeCuts;
  TString AntiCascadeCutsName = Form("%sAntiCascadeCuts", addon.Data());
  coutputAntiCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputAntiCascadeCuts);


  AliAnalysisDataContainer *coutputCascadeOmegaCuts;
  TString CascadeOmegaCutsName = Form("%sCascadeOmegaCuts", addon.Data());
  coutputCascadeOmegaCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeOmegaCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeOmegaCutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputCascadeOmegaCuts);

  AliAnalysisDataContainer *coutputAntiCascadeOmegaCuts;
  TString AntiCascadeOmegaCutsName = Form("%sAntiCascadeOmegaCuts", addon.Data());
  coutputAntiCascadeOmegaCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeOmegaCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeOmegaCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiCascadeOmegaCuts);



  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults", addon.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 9, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("%sResultQA", addon.Data());
  coutputResultQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultQAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultQA);

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample", addon.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultsSample);

  AliAnalysisDataContainer *coutputResultQASample;
  TString ResultQASampleName = Form("%sResultQASample", addon.Data());
  coutputResultQASample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultQASampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQASampleName.Data()));
  mgr->ConnectOutput(task, 12, coutputResultQASample);

  //Cascade Tree:
  AliAnalysisDataContainer *coutputTreeCascade;
  TString TreeCascadeName = Form("%sTreeCascade",addon.Data());
  coutputTreeCascade = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    TreeCascadeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), TreeCascadeName.Data()));
  mgr->ConnectOutput(task, 13, coutputTreeCascade);

  //omega tree:
  AliAnalysisDataContainer *coutputTreeOmega;
  TString TreeOmegaName = Form("%sTreeOmega",addon.Data());
  coutputTreeOmega = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    TreeOmegaName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), TreeOmegaName.Data()));
  mgr->ConnectOutput(task, 14, coutputTreeOmega);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC", addon.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC", addon.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiTrkCutsMC);

//skip v0's as well in MC
/*
    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC", addon.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 17, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC", addon.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 18, coutputAntiv0CutsMC);

*/

    AliAnalysisDataContainer *coutputXiCutsMC;
    TString XiCutsMCName = Form("%sXiCutsMC", addon.Data());
    coutputXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        XiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsMCName.Data()));
    mgr->ConnectOutput(task, 19, coutputXiCutsMC);

    AliAnalysisDataContainer *coutputAntiXiCutsMC;
    TString AntiXiCutsMCName = Form("%sAntiXiCutsMC", addon.Data());
    coutputAntiXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiXiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 20, coutputAntiXiCutsMC);

    AliAnalysisDataContainer *coutputXiOmegaCutsMC;
    TString XiOmegaCutsMCName = Form("%sXiOmegaCutsMC", addon.Data());
    coutputXiOmegaCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        XiOmegaCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiOmegaCutsMCName.Data()));
    mgr->ConnectOutput(task, 17, coutputXiOmegaCutsMC);

    AliAnalysisDataContainer *coutputAntiXiOmegaCutsMC;
    TString AntiXiOmegaCutsMCName = Form("%sAntiXiOmegaCutsMC", addon.Data());
    coutputAntiXiOmegaCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiXiOmegaCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiOmegaCutsMCName.Data()));
    mgr->ConnectOutput(task, 18, coutputAntiXiOmegaCutsMC);
   }




  return task;
}

