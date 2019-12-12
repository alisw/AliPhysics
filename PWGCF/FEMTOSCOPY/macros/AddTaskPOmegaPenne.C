#include "AliAnalysisTaskPOmegaPenne.h"

AliAnalysisTaskPOmegaPenne *AddTaskPOmegaPenne(
    bool isMC = false,
    TString CentEst = "kHM")
{
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

    AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
    evtCuts->CleanUpMult(false, false, false, true);


    // ############# Proton Cuts ###############
    //
    // AliFemtoDreamTrackCuts *TrackCutsProton = new AliFemtoDreamTrackCuts();
    AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts( isMC, true, false, false);
    // TrackCutsProton->SetPlotDCADist(false);
    // TrackCutsProton->SetPlotCombSigma(false);
    // TrackCutsProton->SetIsMonteCarlo(isMC);
    TrackCutsProton->SetCutCharge(1);
    TrackCutsProton->SetFilterBit(128);
    // TrackCutsProton->SetPtRange(0.5, 4.05);
    // TrackCutsProton->SetEtaRange(-0.8, 0.8);
    // TrackCutsProton->SetNClsTPC(80);
    // TrackCutsProton->SetDCAReCalculation(true);
    // TrackCutsProton->SetDCAVtxZ(0.2);
    // TrackCutsProton->SetDCAVtxXY(0.1);
    // TrackCutsProton->SetCutSharedCls(true);
    // TrackCutsProton->SetCutTPCCrossedRows(true, 70, 0.83);
    // TrackCutsProton->SetPID(AliPID::kProton, 0.75);
    // TrackCutsProton->SetCutSmallestSig(true);
    

    // ############# hier anti-Protonen Cuts ###############
    //
    AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts( isMC, true, false, false);
    // AliFemtoDreamTrackCuts *TrackCutsAntiProton = new AliFemtoDreamTrackCuts();
    // TrackCutsAntiProton->SetPlotDCADist(false);
    // TrackCutsAntiProton->SetPlotCombSigma(false);
    // TrackCutsAntiProton->SetIsMonteCarlo(isMC);
    TrackCutsAntiProton->SetCutCharge(-1);
    TrackCutsAntiProton->SetFilterBit(128);
    // TrackCutsAntiProton->SetPtRange(0.5, 4.05);
    // TrackCutsAntiProton->SetEtaRange(-0.8, 0.8);
    // TrackCutsAntiProton->SetNClsTPC(80);
    // TrackCutsAntiProton->SetDCAReCalculation(true);
    // TrackCutsAntiProton->SetDCAVtxZ(0.2);
    // TrackCutsAntiProton->SetDCAVtxXY(0.1);
    // TrackCutsAntiProton->SetCutSharedCls(true);
    // TrackCutsAntiProton->SetCutTPCCrossedRows(true, 70, 0.83);
    // TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75);
    // TrackCutsAntiProton->SetCutSmallestSig(true);


    // ############# Kaon Cuts ###############
    //
    AliFemtoDreamTrackCuts *TrackCutsKaon = AliFemtoDreamTrackCuts::PrimKaonCuts( isMC, true, false, false);
    // AliFemtoDreamTrackCuts *TrackCutsKaon = new AliFemtoDreamTrackCuts();
    // TrackCutsKaon->SetPlotDCADist(false);
    // TrackCutsKaon->SetPlotCombSigma(false);
    // TrackCutsKaon->SetIsMonteCarlo(isMC);
    TrackCutsKaon->SetCutCharge(1);
    TrackCutsKaon->SetFilterBit(128);
    // TrackCutsKaon->SetPtRange(0.15, 999);
    // TrackCutsKaon->SetEtaRange(-0.8, 0.8);
    // TrackCutsKaon->SetNClsTPC(80);
    // TrackCutsKaon->SetDCAReCalculation(true);
    // TrackCutsKaon->SetDCAVtxZ(0.2);
    // TrackCutsKaon->SetDCAVtxXY(0.1);
    // TrackCutsKaon->SetCutSharedCls(true);
    // TrackCutsKaon->SetCutTPCCrossedRows(true, 70, 0.80);
    // TrackCutsKaon->SetPID(AliPID::kKaon, 0.4, 5);
    // TrackCutsKaon->SetCutSmallestSig(true);


    // ############# AntiKaon Cuts ###############
    //
    AliFemtoDreamTrackCuts *TrackCutsAntiKaon = AliFemtoDreamTrackCuts::PrimKaonCuts( isMC, true, false, false);
    // AliFemtoDreamTrackCuts *TrackCutsAntiKaon = new AliFemtoDreamTrackCuts();
    // TrackCutsAntiKaon->SetPlotDCADist(false);
    // TrackCutsAntiKaon->SetPlotCombSigma(false);
    // TrackCutsAntiKaon->SetIsMonteCarlo(isMC);
    TrackCutsAntiKaon->SetCutCharge(-1);
    TrackCutsAntiKaon->SetFilterBit(128);
    // TrackCutsAntiKaon->SetPtRange(0.15, 999);
    // TrackCutsAntiKaon->SetEtaRange(-0.8, 0.8);
    // TrackCutsAntiKaon->SetNClsTPC(80);
    // TrackCutsAntiKaon->SetDCAReCalculation(true);
    // TrackCutsAntiKaon->SetDCAVtxZ(0.2);
    // TrackCutsAntiKaon->SetDCAVtxXY(0.1);
    // TrackCutsAntiKaon->SetCutSharedCls(true);
    // TrackCutsAntiKaon->SetCutTPCCrossedRows(true, 70, 0.80);
    // TrackCutsAntiKaon->SetPID(AliPID::kKaon, 0.4, 5);
    // TrackCutsAntiKaon->SetCutSmallestSig(true);

    std::vector<int> PDGParticles;
    PDGParticles.push_back(2212);   // Protons
    PDGParticles.push_back(2212);

    PDGParticles.push_back(321);    // Kaons
    PDGParticles.push_back(321);

    // vector( size_type count, const T& value, const Allocator& alloc = Allocator());
    std::vector<int> NBins = std::vector<int>(10, 750);
    std::vector<float> kMin = std::vector<float>(10, 0.);
    std::vector<float> kMax = std::vector<float>(10, 3.);
    std::vector<int> pairQA = std::vector<int>(10, 0);
    std::vector<bool> closeRejection = std::vector<bool>(10, false);
    
    pairQA[0] = 11; //p-p
    pairQA[4] = 11; //ap-ap
    closeRejection[0] = true;
    closeRejection[4] = true;

    pairQA[2] = 13;     // ??? p-K
    pairQA[6] = 13;     // ??? ap-ak 

    // ZVtx bins
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
    // Multiplicity bins
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

    AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");
    config->SetZBins(ZVtxBins);
    config->SetMultBins(MultBins);
    config->SetMultBinning(true);
    config->SetPDGCodes(PDGParticles);
    config->SetNBinsHist(NBins);
    config->SetMinKRel(kMin);
    config->SetMaxKRel(kMax);
    config->SetMixingDepth(10);
    config->SetExtendedQAPairs(pairQA);
    config->SetClosePairRejection(closeRejection);
    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);
    
    // task creation
    AliAnalysisTaskPOmegaPenne *task = new AliAnalysisTaskPOmegaPenne("FemtoDreamDefault", isMC);
    if (CentEst == "kInt7")
    {
        task->SelectCollisionCandidates(AliVEvent::kINT7);
        std::cout << "Added kINT7 Trigger \n";
    }
    else if (CentEst == "kHM")
    {
        task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
        std::cout << "Added kHighMultV0 Trigger \n";
    }
    else
    {
        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl;
        std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl;
    }

    task->SetEventCuts(evtCuts);
    task->SetTrackCutsProton(TrackCutsProton);
    task->SetTrackCutsAntiProton(TrackCutsAntiProton);
    task->SetTrackCutsKaon(TrackCutsKaon);
    task->SetTrackCutsAntiKaon(TrackCutsAntiKaon);
    task->SetCollectionConfig(config);

    mgr->AddTask(task);

    TString file = AliAnalysisManager::GetCommonFileName();

    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

    mgr->ConnectInput(task, 0, cinput);

    AliAnalysisDataContainer *coutputQA;
    TString QAName = Form("MyTask");
    coutputQA = mgr->CreateContainer(
        QAName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), QAName.Data()));
    mgr->ConnectOutput(task, 1, coutputQA);

    return task;
}
