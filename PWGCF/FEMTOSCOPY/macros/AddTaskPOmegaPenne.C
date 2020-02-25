#include "AliAnalysisTaskPOmegaPenne.h"
#include <vector>
#include "AliFemtoDreamCascadeCuts.h"

AliAnalysisTaskPOmegaPenne *AddTaskPOmegaPenne( bool isMC = false, TString CentEst = "kHM")
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

    //  ########################### CUTS ##############################
    //
    // Event Cuts
    AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
    evtCuts->CleanUpMult(false, false, false, true);
    evtCuts->SetMultVsCentPlots(true);
    // Track Cuts
    AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    TrackCutsProton->SetFilterBit(128);
    TrackCutsProton->SetCutCharge(1);

    AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    TrackCutsAntiProton->SetFilterBit(128);
    TrackCutsAntiProton->SetCutCharge(-1);

    //  #### Cascade Cuts
    //
    // Xion Cascade
    //
    AliFemtoDreamCascadeCuts *CascadeCutsXion = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
    CascadeCutsXion->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    XiNegCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering
    AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    XiPosCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering
    AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    XiBachCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering

    CascadeCutsXion->Setv0Negcuts(XiNegCuts);
    CascadeCutsXion->Setv0PosCuts(XiPosCuts);
    CascadeCutsXion->SetBachCuts(XiBachCuts);
    CascadeCutsXion->SetPDGCodeCasc(3312);
    CascadeCutsXion->SetPDGCodev0(3122);
    CascadeCutsXion->SetPDGCodePosDaug(2212);
    CascadeCutsXion->SetPDGCodeNegDaug(-211);
    CascadeCutsXion->SetPDGCodeBach(-211);

    // Anti Xion Cascade
    //
    AliFemtoDreamCascadeCuts *AntiCascadeCutsXion = AliFemtoDreamCascadeCuts::XiCuts(isMC, false); 
    AntiCascadeCutsXion->SetXiCharge(1);

    AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    AntiXiNegCuts->SetCutCharge(-1);
    
    AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    AntiXiPosCuts->SetCutCharge(1);
    
    AliFemtoDreamTrackCuts *AntiXiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    AntiXiBachCuts->SetCutCharge(1);

    AntiCascadeCutsXion->Setv0Negcuts(AntiXiNegCuts);
    AntiCascadeCutsXion->Setv0PosCuts(AntiXiPosCuts);
    AntiCascadeCutsXion->SetBachCuts(AntiXiBachCuts);
    AntiCascadeCutsXion->SetPDGCodeCasc(-3312);
    AntiCascadeCutsXion->SetPDGCodev0(-3122);
    AntiCascadeCutsXion->SetPDGCodePosDaug(211);
    AntiCascadeCutsXion->SetPDGCodeNegDaug(-2212);
    AntiCascadeCutsXion->SetPDGCodeBach(211);

    // std::cout << "Hier ist noch /'alles/' OK. Nr.: " << iOkCounter++ << std::endl;
    // evtCuts->SetMinimalBooking(true);
    // TrackCutsProton->SetMinimalBooking(true);
    // TrackCutsAntiProton->SetMinimalBooking(true);
    // CascadeCutsXion->SetMinimalBooking(true);
    // AntiCascadeCutsXion->SetMinimalBooking(true);

    std::vector<int> PDGParticles;
    PDGParticles.push_back(2212);   // Protons
    PDGParticles.push_back(2212);

    PDGParticles.push_back(3312);   // Xions
    PDGParticles.push_back(3312);

    /* std::vector( size_type count, const T& value, const Allocator& alloc = Allocator()); */
    std::vector<int> NBins = std::vector<int>(10, 750);
    std::vector<float> kMin = std::vector<float>(10, 0.);
    std::vector<float> kMax = std::vector<float>(10, 3.);
    std::vector<int> pairQA = std::vector<int>(10, 0);
    std::vector<bool> closeRejection = std::vector<bool>(10, false);
    
    //pairs:
    //pp                0
    //p bar p           1
    //p Xi              2
    //p bar Xi          3
    //bar p bar p       4
    //bar p Xi          5
    //bar p bar Xi      6
    //Xi Xi             7
    //Xi bar Xi         8
    //bar Xi bar Xi     9
    pairQA[0] = 11; 
    pairQA[4] = 11; 
    closeRejection[0] = true;
    closeRejection[4] = true;

    pairQA[2] = 13;     
    pairQA[6] = 13;     
    
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
    MultBins.push_back(84);
    MultBins.push_back(88);
    MultBins.push_back(92);
    MultBins.push_back(96);
    MultBins.push_back(100);

    AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");
    config->SetZBins(ZVtxBins);
    config->SetMultBins(MultBins);
    config->SetMultBinning(true);
    config->SetmTBinning(true);

    config->SetdPhidEtaPlotsSmallK(false);
    config->SetdPhidEtaPlots(false);
    config->SetPhiEtaBinnign(false);
    
    config->SetPDGCodes(PDGParticles);
    config->SetNBinsHist(NBins);
    config->SetMinKRel(kMin);
    config->SetMaxKRel(kMax);
    config->SetMixingDepth(10);
    config->SetExtendedQAPairs(pairQA);
    config->SetClosePairRejection(closeRejection);
    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    // full blast QA
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMassQA(true);

    // task creation
    AliAnalysisTaskPOmegaPenne *task = new AliAnalysisTaskPOmegaPenne("FemtoDreamPOmegaPenne", isMC);
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
    task->SetTrackCutsXion(CascadeCutsXion);
    task->SetTrackCutsAntiXion(AntiCascadeCutsXion);
    task->SetCollectionConfig(config);

    mgr->AddTask(task);

    TString file = AliAnalysisManager::GetCommonFileName();

    // AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    //
    // mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    AliAnalysisDataContainer *coutputQA;
    TString QAName = Form("MyTask");
    coutputQA = mgr->CreateContainer(
                        QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), QAName.Data()));
    mgr->ConnectOutput(task, 1, coutputQA);

    return task;
}
