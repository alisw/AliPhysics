#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE* AddTaskFemtoGranma(
    bool isMC = false,//1
    bool Systematic = false, //2
    TString CentEst = "kInt7",//3
    bool DCAPlots = false,//4
    bool CPAPlots = false,//5
    bool MomReso = false,//6 to set to true only when running on MC
    bool etaPhiPlotsAtTPCRadii=true,//7 to set to true only when running on MC but very Mem. Consuming
    bool CombSigma = false,//8
    bool PileUpRej=true,//9
    bool dPhidEtaPlots=true,//10
    bool ContributionSplitting = false,//11
    bool InvMassPairs=false, //12
    bool kTCentBins=false,//13
    bool DeltaEtaDeltaPhiCut=false,//14
    bool DoSphericityCuts=false, //15
    bool DoSpherocityCuts=false, //16
    const char *swuffix = "",//17
    const char *swuffixvar = "0") {


      // 1    2     3     4     5     6     7    8    9      10   11     12   13    14    15    16   17
      //true,true,false,false,false,false,false,true,false,false,true,false,true,false,false,false,true

  TString suffix=Form("%s",swuffix);
  TString suffixvar=Form("%s",swuffixvar);
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

  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
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
  evtCuts->SetMultVsCentPlots(true);
  evtCuts->SetDoSphericityCuts(DoSphericityCuts);
  evtCuts->SetDoSpherocityCuts(DoSpherocityCuts);
  if(isMC && CentEst=="kHM"){
    evtCuts->SetMultiplicityPercentileMax(5);
  }

  if(DoSphericityCuts){
  if (suffix=="1") {
    evtCuts->SetSphericityCuts(0.,0.3);
  }
  if (suffix=="2") {
    evtCuts->SetSphericityCuts(0.3,0.7);
  }
  if (suffix=="3") {
    evtCuts->SetSphericityCuts(0.7,1.0);
  }
  //4 is reserved to the full sphericity range
  if (suffix=="4") {
    evtCuts->SetSphericityCuts(0.,1.0);
  }
  if (suffix=="5") {
    evtCuts->SetSphericityCuts(0.8,1.0);
  }
  if (suffix=="6") {
    evtCuts->SetSphericityCuts(0.9,1.0);
  }
}else if(!DoSphericityCuts && !DoSpherocityCuts){
    suffix="8";
}

if(DoSpherocityCuts){
if (suffix=="01") {
  evtCuts->SetSpherocityCuts(0.,0.3);
}
if (suffix=="02") {
  evtCuts->SetSpherocityCuts(0.3,0.7);
}
if (suffix=="03") {
  evtCuts->SetSpherocityCuts(0.7,1.0);
}
//4 is reserved to the full sphericity range
if (suffix=="04") {
  evtCuts->SetSpherocityCuts(0.,1.0);
}
if (suffix=="05") {
  evtCuts->SetSpherocityCuts(0.8,1.0);
}
if (suffix=="06") {
  evtCuts->SetSpherocityCuts(0.9,1.0);
}
}else if(!DoSpherocityCuts && !DoSphericityCuts){
  suffix="08";
}

  AliAnalysisTaskGrandma *task = new AliAnalysisTaskGrandma("myFirstTask",
                                                            isMC);

//Track cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, DCAPlots, CombSigma, ContributionSplitting);
  TrackCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, DCAPlots, CombSigma, ContributionSplitting);
  AntiTrackCuts->SetCutCharge(-1);
  //V0 cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(
      isMC,CPAPlots,ContributionSplitting);
  AliFemtoDreamTrackCuts *Posv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
          isMC,PileUpRej,false);
  AliFemtoDreamTrackCuts *Negv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
          isMC,PileUpRej,false);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);//Proton
  v0Cuts->SetPDGCodeNegDaug(211);//Pion
  v0Cuts->SetPDGCodev0(3122);//Lambda
  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(
      isMC, CPAPlots,ContributionSplitting);
  AliFemtoDreamTrackCuts *PosAntiv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
          isMC,PileUpRej,false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
          isMC,PileUpRej,false);
  NegAntiv0Daug->SetCutCharge(-1);
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);//Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);//Proton
  Antiv0Cuts->SetPDGCodev0(-3122);//Lambda


  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");

  std::vector<int> PairQA;
  PairQA.push_back(0);        // p p
  PairQA.push_back(11);         // p barp
  PairQA.push_back(0);        // p Lambda
  PairQA.push_back(12);         // p barLambda
  // PairQA.push_back(0);         // p Xi
  // PairQA.push_back(0);         // p barXi
  PairQA.push_back(0);        // barp barp
  PairQA.push_back(12);         // barp Lambda
  PairQA.push_back(0);        // barp barLambda
  // PairQA.push_back(0);         // barp Xi
  // PairQA.push_back(0);         // barp barXi
  PairQA.push_back(0);         // Lambda Lambda
  PairQA.push_back(22);         // Lambda barLambda
  // PairQA.push_back(0);         // Lambda Xi
  // PairQA.push_back(0);         // Lambda barXi
  PairQA.push_back(0);         // barLambda barLamb
  // PairQA.push_back(0);         // barLambda Xi
  // PairQA.push_back(0);         // barLambda barXi
  // PairQA.push_back(0);         // Xi Xi
  // PairQA.push_back(0);         // Xi barXi
  // PairQA.push_back(0);         // barXi barXi
  config->SetExtendedQAPairs(PairQA);

  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);
  config->SetPDGCodes(PDGParticles);

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
  config->SetZBins(ZVtxBins);

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

  std::vector<int> centBins;
  centBins.push_back(20);
  centBins.push_back(40);
  centBins.push_back(90);
  config->SetCentBins(centBins);
  config->SetkTCentralityBinning(kTCentBins);


if(isMC)
  {
  config->SetdPhidEtaPlots(dPhidEtaPlots);  // warsaw like plots
  std::cout<<"in MC to make dETAdPHI plots"<<dPhidEtaPlots<<std::endl;
} else{
  config->SetdPhidEtaPlots(dPhidEtaPlots);  // warsaw like plots
}

if (MomReso) {
  if (isMC) {
    config->SetMomentumResolution(true);//kstar true vs. kstar reco
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


  std::vector<int> NBins;
  NBins.push_back(750);  // p p
  NBins.push_back(750);  // p barp
  NBins.push_back(750);  // p Lambda
  NBins.push_back(750);  // p barLambda
//  NBins.push_back(750);  // p Xi
//  NBins.push_back(750);  // p barXi
  NBins.push_back(750);  // barp barp
  NBins.push_back(750);  // barp Lambda
  NBins.push_back(750);  // barp barLambda
//  NBins.push_back(750);  // barp Xi
//  NBins.push_back(750);  // barp barXi
  NBins.push_back(750);  // Lambda Lambda
  NBins.push_back(750);  // Lambda barLambda
//  NBins.push_back(750);  // Lambda Xi
//  NBins.push_back(750);  // Lambda barXi
  NBins.push_back(750);  // barLambda barLambda
//  NBins.push_back(750);  // barLambda Xi
//  NBins.push_back(750);  // barLambda barXi
//  NBins.push_back(750);  // Xi Xi
//  NBins.push_back(750);  // Xi barXi
//  NBins.push_back(750);  // barXi barXi
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
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
//  kMin.push_back(0.);
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
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
//  kMax.push_back(3.);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMultBinning(true);
  config->SetUseEventMixing(true);
  config->SetMixingDepth(10);
  config->SetCentBinning(false);
  config->SetkTBinning(false);
  config->SetmTBinning(true);


  config->SetUsePhiSpinning(false);
//  config->SetSpinningDepth(10);
//  config->SetUseStravinskyMethod(false);

  config->SetMinimalBookingME(false);
  config->SetMinimalBookingSample(true);

  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);


  if (CentEst == "kInt7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kMB") {
    task->SelectCollisionCandidates(AliVEvent::kMB);
    std::cout << "Added kMB Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
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

  if (Systematic) {
    if(suffixvar == "0"){

      printf("Running Default Cut Variations\n");
    }
   if (suffixvar == "1") {
     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "2") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

   } else if (suffixvar == "3") {
     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "4") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "5") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "6") {
     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

   } else if (suffixvar == "7") {
     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "8") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

   } else if (suffixvar == "9") {
     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "10") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "11") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "12") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

   } else if (suffixvar == "13") {
     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "14") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

   } else if (suffixvar == "15") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "16") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "17") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "18") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

   } else if (suffixvar == "19") {
     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "20") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

   } else if (suffixvar == "21") {

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

   } else if (suffixvar == "22") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

   } else if (suffixvar == "23") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "24") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     //XI

   } else if (suffixvar == "25") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.83, 0.83);
     AntiTrackCuts->SetEtaRange(-0.83, 0.83);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "26") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

   } else if (suffixvar == "27") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "28") {
     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

   } else if (suffixvar == "29") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
   } else if (suffixvar == "30") {
     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "31") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "32") {

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
   } else if (suffixvar == "33") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
   } else if (suffixvar == "34") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

   } else if (suffixvar == "35") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "36") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "37") {
     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

   } else if (suffixvar == "38") {
     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
   } else if (suffixvar == "39") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     TrackCuts->SetNClsTPC(90);
     AntiTrackCuts->SetNClsTPC(90);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "40") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetNClsTPC(70);
     AntiTrackCuts->SetNClsTPC(70);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "41") {
     TrackCuts->SetPtRange(0.4, 4.05);
     AntiTrackCuts->SetPtRange(0.4, 4.05);

     TrackCuts->SetEtaRange(-0.77, 0.77);
     AntiTrackCuts->SetEtaRange(-0.77, 0.77);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetEtaRange(-0.77, 0.77);
     Negv0Daug->SetEtaRange(-0.77, 0.77);
     PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
     NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "42") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     v0Cuts->SetCutDCADaugTov0Vtx(1.2);
     Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
   } else if (suffixvar == "43") {
     TrackCuts->SetPtRange(0.6, 4.05);
     AntiTrackCuts->SetPtRange(0.6, 4.05);

     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
   } else if (suffixvar == "44") {
     TrackCuts->SetEtaRange(-0.85, 0.85);
     AntiTrackCuts->SetEtaRange(-0.85, 0.85);

     TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
     AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

     config->SetDeltaEtaMax(0.012);
     config->SetDeltaPhiMax(0.012);

     v0Cuts->SetCutCPA(0.995);
     Antiv0Cuts->SetCutCPA(0.995);

     Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
     Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
     NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

     Posv0Daug->SetNClsTPC(80);
     Negv0Daug->SetNClsTPC(80);
     PosAntiv0Daug->SetNClsTPC(80);
     NegAntiv0Daug->SetNClsTPC(80);

     Posv0Daug->SetEtaRange(-0.83, 0.83);
     Negv0Daug->SetEtaRange(-0.83, 0.83);
     PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
     NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

     v0Cuts->SetCutDCADaugToPrimVtx(0.06);
     Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
   }
 }


  task->SetEvtCutQA(true);
  task->SetTrackBufferSize(2000);
  task->SetEventCuts(evtCuts);
  task->SetTrackCuts(TrackCuts);
  task->SetAntiTrackCuts(AntiTrackCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);
  TString addon="";
  if (CentEst=="kInt7") {
    addon+="MB";
  } else if (CentEst=="kHM") {
    addon+="HM";
  }

  AliAnalysisDataContainer *coutputQA;

  std::cout << "CONTAINTER NAME: " << addon.Data() << " " << suffix.Data() << std::endl;
  TString QAName;
  if(Systematic){
     QAName = Form("%sQA%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
  }else{
   QAName = Form("%sQA%s",addon.Data(),suffix.Data());
  }
  coutputQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      QAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  AliAnalysisDataContainer *coutputEvtCuts;
  TString EvtCutsName;
  if(Systematic){
   EvtCutsName = Form("%sEvtCuts%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
 }else{
   EvtCutsName = Form("%sEvtCuts%s",addon.Data(),suffix.Data());
 }
  coutputEvtCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      EvtCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *couputTrkCuts;
  TString TrackCutsName;
if(Systematic){
  TrackCutsName = Form("%sTrackCuts%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
}else{
 TrackCutsName = Form("%sTrackCuts%s",addon.Data(),suffix.Data());
}
  couputTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      TrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  AliAnalysisDataContainer *coutputAntiTrkCuts;
  TString AntiTrackCutsName;
  if(Systematic){
  AntiTrackCutsName = Form("%sAntiTrackCuts%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
}else{
  AntiTrackCutsName = Form("%sAntiTrackCuts%s",addon.Data(),suffix.Data());
}
  coutputAntiTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiTrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  AliAnalysisDataContainer *couputv0Cuts;
  TString v0CutsName;
  if(Systematic){
    v0CutsName = Form("%sv0Cuts%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
  } else{
   v0CutsName = Form("%sv0Cuts%s",addon.Data(),suffix.Data());
}
  couputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, couputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName;
if(Systematic){
  Antiv0CutsName = Form("%sAntiv0Cuts%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
}else{
   Antiv0CutsName = Form("%sAntiv0Cuts%s",addon.Data(),suffix.Data());
}
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName;
  if(Systematic){
    ResultsName = Form("%sResults%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
  }else{
   ResultsName = Form("%sResults%s",addon.Data(),suffix.Data());
}
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 7, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName;
  if(Systematic){
    ResultQAName = Form("%sResultQA%s_%s",addon.Data(),suffix.Data(),suffixvar.Data());
  }else{
  ResultQAName = Form("%sResultQA%s",addon.Data(),suffix.Data());
}
  coutputResultQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultQAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 8, coutputResultQA);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC%s",addon.Data(),suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 9, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s",addon.Data(),suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s",addon.Data(),suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s",addon.Data(),suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputAntiv0CutsMC);
  }

  return task;

}
