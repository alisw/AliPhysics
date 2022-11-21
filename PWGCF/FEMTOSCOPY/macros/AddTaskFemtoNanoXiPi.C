#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoXiPi.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskFemtoNanoXiPi(bool isMC = false,              // 2
                                        int fFilterBit = 128,           // 3
                                        bool Systematic = false,        // 6
                                        bool DoAncestors = false,       // 10
                                        const char *cutVariation = "0") // 11
{
    TString suffix = TString::Format("%s", cutVariation);

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        Error("AddTaskFemtoNanoXiPi()", "No analysis manager found.");
        return 0x0;
    }

    // ================== GetInputEventHandler =============================
    AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

    //========= Init subtasks and start analysis ============================
    // Event Cuts
    AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
    //  void CleanUpMult(bool SPD, bool v0A, bool v0C, bool RefMult) {
    evtCuts->CleanUpMult(false, false, false, true);

    /// Pion cuts
    AliFemtoDreamTrackCuts *TrackCutsPion = NULL;
    AliFemtoDreamTrackCuts *TrackCutsAntiPion = NULL;

    TrackCutsPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
    TrackCutsPion->SetFilterBit(fFilterBit);
    TrackCutsPion->SetCutCharge(1);
    TrackCutsAntiPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
    TrackCutsAntiPion->SetFilterBit(fFilterBit);
    TrackCutsAntiPion->SetCutCharge(-1);

    // Cascade Cuts
    AliFemtoDreamCascadeCuts *CascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
        isMC, false);
    CascadeCuts->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
        isMC, true, false);
    XiNegCuts->SetCheckTPCRefit(false); // for nanos this is already done while prefiltering
    AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
        isMC, true, false);
    XiPosCuts->SetCheckTPCRefit(false); // for nanos this is already done while prefiltering
    AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(
        isMC, true, false);
    XiBachCuts->SetCheckTPCRefit(false); // for nanos this is already done while prefiltering

    CascadeCuts->Setv0Negcuts(XiNegCuts);
    CascadeCuts->Setv0PosCuts(XiPosCuts);
    CascadeCuts->SetBachCuts(XiBachCuts);
    CascadeCuts->SetPDGCodeCasc(3312);
    CascadeCuts->SetPDGCodev0(3122);
    CascadeCuts->SetPDGCodePosDaug(2212);
    CascadeCuts->SetPDGCodeNegDaug(-211);
    CascadeCuts->SetPDGCodeBach(-211);

    AliFemtoDreamCascadeCuts *AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
        isMC, false);
    AntiCascadeCuts->SetXiCharge(1);
    AliFemtoDreamTrackCuts *AntiXiNegCuts =
        AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    AntiXiNegCuts->SetCutCharge(-1);
    AntiXiNegCuts->SetCheckTPCRefit(false); // for nanos this is already done while prefiltering
    AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
        isMC, true, false);
    AntiXiPosCuts->SetCutCharge(1);
    AntiXiPosCuts->SetCheckTPCRefit(false); // for nanos this is already done while prefiltering
    AliFemtoDreamTrackCuts *AntiXiBachCuts =
        AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    AntiXiBachCuts->SetCutCharge(1);
    AntiXiBachCuts->SetCheckTPCRefit(false); // for nanos this is already done while prefiltering

    AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
    AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
    AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
    AntiCascadeCuts->SetPDGCodeCasc(-3312);
    AntiCascadeCuts->SetPDGCodev0(-3122);
    AntiCascadeCuts->SetPDGCodePosDaug(211);
    AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
    AntiCascadeCuts->SetPDGCodeBach(211);

    if (Systematic)
    {
        evtCuts->SetMinimalBooking(true);
        TrackCutsPion->SetMinimalBooking(true);
        TrackCutsAntiPion->SetMinimalBooking(true);
        CascadeCuts->SetMinimalBooking(true);
        AntiCascadeCuts->SetMinimalBooking(true);
    }

    AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                  "Femto", false);

    // Femto Collection
    std::vector<int> PDGParticles;
    PDGParticles.push_back(211);  // pion +
    PDGParticles.push_back(-211); // pion -
    PDGParticles.push_back(3312); // Cascade
    PDGParticles.push_back(-3312);

    std::vector<int> NBins;
    std::vector<float> kMin;
    std::vector<float> kMax;
    std::vector<int> pairQA;
    std::vector<int> pairQASyst;
    std::vector<bool> closeRejection;
    // pairs:
    // pi+ pi+           0
    // pi+ pi-           1

    // pi+ Xi-           2 <--------
    // pi+ AntiXi+       3

    // pi- pi-           4
    // pi- Xi-           5
    // pi- AntiXi+       6 <---------

    // Xi Xi             7
    // Xi AntiXi         8
    // AntiXi AntiXi     9

    const int nPairs = 10;
    for (int i = 0; i < nPairs; ++i)
    {
        pairQA.push_back(0);
        closeRejection.push_back(false);
        NBins.push_back(1500);
        kMin.push_back(0.);
        kMax.push_back(6.);
    }

    if (Systematic)
    {
        pairQA[2] = 12;
        pairQA[6] = 12;
        closeRejection[2] = true; // pp
        closeRejection[6] = true; // barp barp
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

        // pairQA[7] = 22;
        // pairQA[8] = 22;
        // pairQA[9] = 22;
        closeRejection[0] = true;
        closeRejection[1] = true;
        closeRejection[2] = true;
        closeRejection[3] = true;
        closeRejection[4] = true;
        closeRejection[5] = true;
        closeRejection[6] = true;
        // closeRejection[7] = true;
        // closeRejection[8] = true;
        // closeRejection[9] = true;
    }

    config->SetPDGCodes(PDGParticles);
    config->SetNBinsHist(NBins);
    config->SetMinKRel(kMin);
    config->SetMaxKRel(kMax);
    config->SetClosePairRejection(closeRejection);
    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);
    config->SetExtendedQAPairs(pairQA);

    config->SetMixingDepth(20);
    config->SetUseEventMixing(true);

    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

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

    if (Systematic)
    {
        config->SetMinimalBookingME(true);
    }
    else if (!Systematic)
    {
        config->SetPtQA(true);
        config->SetMassQA(true);
    }

    if (isMC)
    {
        config->SetMomentumResolution(true); // kstar true vs. kstar reco
    }
    else
    {
        std::cout
            << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
    }

    // Common/Non Common Ancestors
    if (isMC && DoAncestors)
    {
        config->SetAncestors(true);
        config->GetDoAncestorsPlots();
    }

    AliAnalysisTaskNanoXiPi *task = new AliAnalysisTaskNanoXiPi("femtoXiPi", isMC);

    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);

    task->SetEventCuts(evtCuts);
    task->SetPionCuts(TrackCutsPion);
    task->SetAntiPionCuts(TrackCutsAntiPion);
    task->SetXiCuts(CascadeCuts);
    task->SetAntiXiCuts(AntiCascadeCuts);
    task->SetCorrelationConfig(config);
    mgr->AddTask(task);

    TString addon = "HM";

    TString file = AliAnalysisManager::GetCommonFileName();

    mgr->ConnectInput(task, 0, cinput);

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

    TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
    AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));
    mgr->ConnectOutput(task, 3, couputTrkCuts);

    TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                     suffix.Data());
    AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
    mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

    AliAnalysisDataContainer *coutputXiCuts;
    TString XiCutsName = Form("%sXiCuts%s", addon.Data(), suffix.Data());
    coutputXiCuts = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        XiCutsName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsName.Data()));
    mgr->ConnectOutput(task, 5, coutputXiCuts);

    AliAnalysisDataContainer *coutputAntiXiCuts;
    TString AntiXiCutsName = Form("%sAntiXiCuts%s", addon.Data(), suffix.Data());
    coutputAntiXiCuts = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiXiCutsName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsName.Data()));
    mgr->ConnectOutput(task, 6, coutputAntiXiCuts);

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

    AliAnalysisDataContainer *coutputResultsSample;
    TString ResultsSampleName = Form("%sResultsSample%s", addon.Data(),
                                     suffix.Data());
    coutputResultsSample = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        ResultsSampleName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), ResultsSampleName.Data()));
    mgr->ConnectOutput(task, 9, coutputResultsSample);

    AliAnalysisDataContainer *coutputResultsSampleQA;
    TString ResultsSampleQAName = Form("%sResultsSampleQA%s", addon.Data(),
                                       suffix.Data());
    coutputResultsSampleQA = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        ResultsSampleQAName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));
    mgr->ConnectOutput(task, 10, coutputResultsSampleQA);

    if (isMC)
    {
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

        AliAnalysisDataContainer *coutputXiCutsMC;
        TString XiCutsMCName = Form("%sXiCutsMC%s", addon.Data(), suffix.Data());
        coutputXiCutsMC = mgr->CreateContainer(
            //@suppress("Invalid arguments") it works ffs
            XiCutsMCName.Data(),
            TList::Class(),
            AliAnalysisManager::kOutputContainer,
            Form("%s:%s", file.Data(), XiCutsMCName.Data()));
        mgr->ConnectOutput(task, 13, coutputXiCutsMC);

        AliAnalysisDataContainer *coutputAntiXiCutsMC;
        TString AntiXiCutsMCName = Form("%sAntiXiCutsMC%s", addon.Data(), suffix.Data());
        coutputAntiXiCutsMC = mgr->CreateContainer(
            //@suppress("Invalid arguments") it works ffs
            AntiXiCutsMCName.Data(),
            TList::Class(),
            AliAnalysisManager::kOutputContainer,
            Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
        mgr->ConnectOutput(task, 14, coutputAntiXiCutsMC);
    }

    return task;
}
