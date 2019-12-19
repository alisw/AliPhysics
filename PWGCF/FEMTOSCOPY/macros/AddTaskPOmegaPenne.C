#include "AliAnalysisTaskPOmegaPenne.h"
#include <vector>
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoXioton.h"
#include "AliAnalysisTaskAODXioton.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

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

    AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
    evtCuts->CleanUpMult(false, false, false, true);


    //  ########################### CUTS ##############################
    //
    //  Proton Cuts
    //
    AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts( isMC, true, false, false);
    TrackCutsProton->SetCutCharge(1);
    // filterbit already set in method to 128

    // anti-Proton Cuts
    //
    AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts( isMC, true, false, false);
    TrackCutsAntiProton->SetCutCharge(-1);
    TrackCutsAntiProton->SetFilterBit(128);

    // Kaon Cuts
    //
    AliFemtoDreamTrackCuts *TrackCutsKaon = AliFemtoDreamTrackCuts::PrimKaonCuts( isMC, true, false, false);
    TrackCutsKaon->SetCutCharge(1);
    TrackCutsKaon->SetFilterBit(128);

    // AntiKaon Cuts
    //
    AliFemtoDreamTrackCuts *TrackCutsAntiKaon = AliFemtoDreamTrackCuts::PrimKaonCuts( isMC, true, false, false);
    TrackCutsAntiKaon->SetCutCharge(-1);
    TrackCutsAntiKaon->SetFilterBit(128);

    // //Cascade Cuts
    // AliFemtoDreamCascadeCuts *CascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
    // CascadeCuts->SetXiCharge(-1);
    // AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    // XiNegCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering
    // AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    // XiPosCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering
    // AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    // XiBachCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering

    // CascadeCuts->Setv0Negcuts(XiNegCuts);
    // CascadeCuts->Setv0PosCuts(XiPosCuts);
    // CascadeCuts->SetBachCuts(XiBachCuts);
    // CascadeCuts->SetPDGCodeCasc(3312);
    // CascadeCuts->SetPDGCodev0(3122);
    // CascadeCuts->SetPDGCodePosDaug(2212);
    // CascadeCuts->SetPDGCodeNegDaug(-211);
    // CascadeCuts->SetPDGCodeBach(-211);

    // AliFemtoDreamCascadeCuts *AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
    // AntiCascadeCuts->SetXiCharge(1);
    // AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    // AntiXiNegCuts->SetCutCharge(-1);
    // AntiXiNegCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering
    // AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    // AntiXiPosCuts->SetCutCharge(1);
    // AntiXiPosCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering
    // AliFemtoDreamTrackCuts *AntiXiBachCuts =
    //     AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    // AntiXiBachCuts->SetCutCharge(1);
    // AntiXiBachCuts->SetCheckTPCRefit(false); //for nanos this is already done while prefiltering

    // AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
    // AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
    // AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
    // AntiCascadeCuts->SetPDGCodeCasc(-3312);
    // AntiCascadeCuts->SetPDGCodev0(-3122);
    // AntiCascadeCuts->SetPDGCodePosDaug(211);
    // AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
    // AntiCascadeCuts->SetPDGCodeBach(211);



    std::vector<int> PDGParticles;
    PDGParticles.push_back(2212);   // Protons
    PDGParticles.push_back(2212);

    PDGParticles.push_back(321);    // Kaons
    PDGParticles.push_back(321);

    // vector( size_type count, const T& value, const Allocator& alloc = Allocator());
    // std::vector<int> NBins = std::vector<int>(10, 750);
    // std::vector<float> kMin = std::vector<float>(10, 0.);
    // std::vector<float> kMax = std::vector<float>(10, 3.);
    // std::vector<int> pairQA = std::vector<int>(10, 0);
    // std::vector<bool> closeRejection = std::vector<bool>(10, false);
    std::vector<int> NBins;
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
    Nbins.push_Back(750);
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
    std::vector<int> pairQA;
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    pairQA.push_back(0);
    std::vector<bool> closeRejection;
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);
    closeRejection.push_back(false);

    pairQA[0] = 11; //p-p
    pairQA[4] = 11; //ap-ap
    closeRejection[0] = true;
    closeRejection[4] = true;

    pairQA[2] = 13;     // p-K
    pairQA[6] = 13;     // ap-ak

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
    task->SetTrackCutsKaon(TrackCutsKaon);
    task->SetTrackCutsAntiKaon(TrackCutsAntiKaon);
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
