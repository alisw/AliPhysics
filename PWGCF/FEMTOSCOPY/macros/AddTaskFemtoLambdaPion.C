#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskLambdaPion.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskFemtoLambdaPion(bool isMC = true,                 // MC run or not
                                              TString CentEst = "kInt7",    // trigger selection, "kInt7" = Minimum bias, "kHighMultV0" high multiplicity triggered by the V0 detector
                                              int filterBit = 128,          // Track selection feature
                                              int WhichPionCut = 0,         // Irrelevant for now
                                              const char *sTcut = "0",      // "0" to avoid spericity cuts, "1" to have them
                                              bool DoAncestors = false,     // only important when running MC
                                              bool IsSystematics = false,   // true to evaluate systematic uncertainties
                                              AliAnalysisTaskLambdaPion::PCSettings pcsettings = AliAnalysisTaskLambdaPion::PCSettings::NoPC,  // choose pair cleaner
                                              bool usenolambdaevt = true,            // true to discard events with neither Lambda or AntiLambda
                                              double dauPIDCut = 2,
                                              const char *cutVariation = "0")
//                                              int binwidth = 1)             // relative bin width for k* histos with respect to 4 MeV/c
{
  TString suffix = TString::Format("%s", cutVariation);
  TString sTsuffix = TString::Format("%s", sTcut);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler())
  {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  //========= Init subtasks and start analysis ============================
  // Event Cuts 
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2(); // uses cuts for events selected by ALICE collaboration
  evtCuts->CleanUpMult(false, false, false, true); // the arguments allow to choose which detectors' multiplicity one wants to consider,
                                                   // first we have SPD, then V0A, then V0C and then Ref08Mult

  if (sTsuffix == "1")
  {
    evtCuts->SetSphericityCuts(0.7,1);
  }

  // Lambda --> p + pi- cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);  // enable PCA plots and Split Contrib, this method sets
                                                                                    // all cuts parameters for pT, charge, DCA, invmass
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false, dauPIDCut); // (bool isMC, bool PileUpRej, bool ContribSplitting)
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false, dauPIDCut);

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);  // saves selected track cuts criteria for charged daughter of neutral particle
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212); // Proton
  v0Cuts->SetPDGCodeNegDaug(211);  // Pion
  v0Cuts->SetPDGCodev0(3122);      // Lambda

  // AntiLambda --> Antip + pi+ cuts
  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false, dauPIDCut);  // Select Pi+ as positive daughter, this means we are considering AntiLambda 
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false, dauPIDCut);
  NegAntiv0Daug->SetCutCharge(-1);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  // Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212); // Proton
  Antiv0Cuts->SetPDGCodev0(-3122);     // Lambda

  // Cuts on correlation pions (Marcel)
  AliFemtoDreamTrackCuts *TrackPosPionCuts = NULL;
  AliFemtoDreamTrackCuts *TrackCutsAntiPion = NULL;

  TrackPosPionCuts = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false); // (bool isMC, bool DCAPlots, bool CombSigma, bool ContribSplitting), sets options for plot
                                                                                     // of DCA distribution, CombSigma, ContribSplitting
  TrackPosPionCuts->SetFilterBit(filterBit);
  TrackPosPionCuts->SetPtRange(0.14, 2.);
  TrackPosPionCuts->SetCutCharge(1);
  TrackCutsAntiPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackCutsAntiPion->SetFilterBit(filterBit);
  TrackCutsAntiPion->SetPtRange(0.14, 2.);
  TrackCutsAntiPion->SetCutCharge(-1);

  if (IsSystematics)
  {
    evtCuts->SetMinimalBooking(true); // minimal booking defines which histograms to save since with systematics many repeat each other
    TrackPosPionCuts->SetMinimalBooking(true);
    TrackCutsAntiPion->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }

  // Correlated particles definition by PDG code
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211);  // 0 Pion Plus
  PDGParticles.push_back(211);  // 1 Pion Minus
  PDGParticles.push_back(3122); // 2 Lambda
  PDGParticles.push_back(3122); // 3 antiLambda

  // We need to set the ZVtx bins
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
  // The Multiplicity bins are set here
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

  // Number of bins of hitograms
  std::vector<int> NBins;
  // minimum k* value
  std::vector<float> kMin;
  // maximum k* value
  std::vector<float> kMax;
  // pair QA extended
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  std::vector<float> mTBins = {0.9, 1.15, 1.4, 4.5, 999.};

  // pairs:
  // Pi+Pi+         0
  // Pi+Pi-         1
  // Pi+ La         2
  // Pi+ bar La     3
  // Pi-Pi-         4
  // Pi- La         5
  // Pi- bar La     6
  // La La          7
  // La La bar      8
  // La bar La bar  9

  const int npairs = 10;
  // settings for k* distribution histos
  for (int i = 0; i < npairs; i++)
  {
    NBins.push_back(1500);
    closeRejection.push_back(false);
    kMin.push_back(0.);
    kMax.push_back(6.);
    pairQA.push_back(0);
  }

  if (IsSystematics)
  {
    pairQA[2] = 12;
    pairQA[3] = 12;
    pairQA[5] = 12;
    pairQA[6] = 12;
  }
  else
  {
    pairQA[0] = 11;
    pairQA[1] = 11;
    pairQA[2] = 12;
    pairQA[3] = 12;
    pairQA[4] = 11;
    pairQA[5] = 12;
    pairQA[6] = 12;
    pairQA[7] = 22;
    pairQA[8] = 22;
    pairQA[9] = 22;
    closeRejection[0] = true;
    closeRejection[1] = true;
    closeRejection[4] = true;
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");

  config->SetPDGCodes(PDGParticles);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);  // save k* distro for each multiplicity bin both for SE and ME
  config->SetClosePairRejection(closeRejection); // true only if we don't calculate systematic effects
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetExtendedQAPairs(pairQA);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetUseEventMixing(true); // ????
  config->SetMixingDepth(30);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  config->SetMultBinning(true);
  config->SetmTBins(mTBins);
  config->SetDomTMultBinning(true);
  config->SetmTBinning(true);

  if (IsSystematics)
  {
    config->SetMinimalBookingME(true);
  }
  else if (!IsSystematics)
  {
    config->SetPtQA(true);
    config->SetMassQA(true);
    config->SetkTBinning(true);
  }

  if (isMC)
  {
    config->SetMomentumResolution(true);
    if (DoAncestors)
    {
      config->SetAncestors(true);
      config->GetDoAncestorsPlots();
    }
  }
  else
  {
    std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

  // Variation cuts
  const float PionPtlow = 0.09;
  const float PionPtup = 0.19;
  const float PionEtaLow = 0.75;
  const float PionEtaUp = 0.85;
  const float PionNClsLow = 70;
  const float PionNClsUp = 90;
  const float PionPtMax = 999;

  AliPID::EParticleType aliPIDParticle;
  aliPIDParticle = AliPID::kPion;
  std::map<std::string, float> PionPIDTight; // tight pion selection criteria
  std::map<std::string, float> PionPIDLoose; // loose pion selection criteria

  PionPIDTight = {
      {"COMB", 2.7},
      {"TPC", 2.7},
      {"EXCLUSION", 3.3},
  }; // for SetPIDkd() when using oton's K selection
  PionPIDLoose = {
      {"COMB", 3.3},
      {"TPC", 3.3},
      {"EXCLUSION", 2.7},
  };

  // Systematic variations (taken from pφ and ΛΚ)
  if (IsSystematics)
  {
    if (suffix == "1")
    {
      TrackPosPionCuts->SetNClsTPC(PionNClsUp);   // Number of TPC clusters
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      v0Cuts->SetCutCPA(0.995);     // CPA stands for cosine of pointing angle
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);     
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);  // DCA cut with respect to the primary vertex
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "2")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);  // (AliPID::EParticleType pid, float pTPCThresh, float sigVal = 3)
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    }
    else if (suffix == "3")
    {
      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "4")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);
      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "5")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      // CHANGE!!!
      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]); 
      // (bool iskaon = true, bool isramona = false, float COMBcut = 3., float TPCcut = 3., float EXCLUSIONcut = 3.)
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "6")
    {
      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
    }
    else if (suffix == "7")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

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
    }
    else if (suffix == "8")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "9")
    {
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
    }
    else if (suffix == "10")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

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
    }
    else if (suffix == "11")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

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
    }
    else if (suffix == "12")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "13")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "14")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "15")
    {
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
    }
    else if (suffix == "16")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
    }
    else if (suffix == "17")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "18")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "19")
    {
      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "20")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "21")
    {

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "22")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    }
    else if (suffix == "23")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

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
    }
    else if (suffix == "24")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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

      // XI
    }
    else if (suffix == "25")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiPion->SetEtaRange(-0.83, 0.83);

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
    }
    else if (suffix == "26")
    {
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "27")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
    }
    else if (suffix == "28")
    {
      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "29")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "30")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

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
    }
    else if (suffix == "31")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "32")
    {

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "33")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

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
    }
    else if (suffix == "34")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "35")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "36")
    {
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
    }
    else if (suffix == "37")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "38")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
    }
    else if (suffix == "39")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
    }
    else if (suffix == "40")
    {
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "41")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "42")
    {
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

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
    }
    else if (suffix == "43")
    {
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "44")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

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

  // Apply sidebands cuts
  if (suffix == "50") {
    // Left sideband 1
    v0Cuts->SetCutWindow(1.080, 1.100);
    Antiv0Cuts->SetCutWindow(1.080, 1.100);
  }
  if (suffix == "51") {
    // Right sideband 1
    v0Cuts->SetCutWindow(1.135, 1.155);
    Antiv0Cuts->SetCutWindow(1.135, 1.155);
  }
  if (suffix == "52") {
    // Right sideband 2
    v0Cuts->SetCutWindow(1.170, 1.200);
    Antiv0Cuts->SetCutWindow(1.170, 1.200);
  }
  if (suffix == "53") {
    // ProtonLambda - Left sideband 1
    v0Cuts->SetCutWindow(1.085, 1.103);
    Antiv0Cuts->SetCutWindow(1.085, 1.103);
  }
  if (suffix == "54") {
    // ProtonLambda - Left sideband 2
    v0Cuts->SetCutWindow(1.090, 1.103);
    Antiv0Cuts->SetCutWindow(1.090, 1.103);
  }
  if (suffix == "55") {
    // ProtonLambda - Left sideband 3
    v0Cuts->SetCutWindow(1.095, 1.108);
    Antiv0Cuts->SetCutWindow(1.095, 1.108);
  }
  if (suffix == "56") {
    // ProtonLambda - Right sideband 1
    v0Cuts->SetCutWindow(1.129, 1.155);
    Antiv0Cuts->SetCutWindow(1.129, 1.155);
  }
  if (suffix == "57") {
    // ProtonLambda - Right sideband 2
    v0Cuts->SetCutWindow(1.129, 1.145);
    Antiv0Cuts->SetCutWindow(1.129, 1.145);
  }
  if (suffix == "58") {
    // ProtonLambda - Right sideband 3
    v0Cuts->SetCutWindow(1.124, 1.140);
    Antiv0Cuts->SetCutWindow(1.124, 1.140);
  }
  if (suffix == "59") {
    // ProtonLambda - Right sideband 4
    v0Cuts->SetCutWindow(1.124, 1.135);
    Antiv0Cuts->SetCutWindow(1.124, 1.135);
  }

  // now we create the task
  AliAnalysisTaskLambdaPion *task =
      new AliAnalysisTaskLambdaPion(
          "AliAnalysisTaskLambdaPion", isMC, pcsettings, usenolambdaevt);

  // trigger selection according to macro arguments
  if (CentEst == "kInt7")
  {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  }
  else if (CentEst == "kHM")
  {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  }
  else
  {
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will "
                 "be empty!"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
  }

  // Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetLambdaCuts(v0Cuts);
  task->SetAntiLambdaCuts(Antiv0Cuts);
  task->SetPosPionCuts(TrackPosPionCuts);
  task->SetNegPionCuts(TrackCutsAntiPion);
  task->SetCollectionConfig(config);

  // mgr adds a user task to the global list of tasks
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  TString addon = "";
  if (CentEst == "kInt7") {
    addon += "MB";
  }
  else if (CentEst == "kHM") {
    addon += "HM";
  }
  // Output file structure definition
  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
      QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sLambdaCuts%s", addon.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiLambdaCuts%s", addon.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiv0Cuts);

  TString TrackCutsName = Form("%sPionPlusCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 5, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sPionMinusCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 7, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 8, coutputResultsQA);

  // other histograms specific for MC
  if (isMC)
  {
    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 9, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputAntiv0CutsMC);

    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputAntiTrkCutsMC);
  }

  return task;
}
