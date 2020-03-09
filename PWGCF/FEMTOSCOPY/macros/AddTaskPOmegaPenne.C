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

    // Event Cuts
    //
    AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
    evtCuts->CleanUpMult(false, false, false, true);
    evtCuts->SetMultVsCentPlots(true);

    // Track Cuts
    //
    AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    TrackCutsProton->SetFilterBit(128);
    TrackCutsProton->SetCutCharge(1);

    AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    TrackCutsAntiProton->SetFilterBit(128);
    TrackCutsAntiProton->SetCutCharge(-1);

    // Lambda Cuts
    //
    AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);

    AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);

    v0Cuts->SetCutTransverseRadius(1.4, 200);           // damit v0s erst ab möglichem xiCut berücksichtigt werden - Xi mittlerer flugweg etwa 4,9 cm
                                                        // dann hau ich auch alle raus d
    v0Cuts->SetCutCPA(0.97);                            // CPA für sekundäre kleiner?
    v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
    v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
    v0Cuts->SetPDGCodePosDaug(2212); //Proton
    v0Cuts->SetPDGCodeNegDaug(211);  //Pion
    v0Cuts->SetPDGCodev0(3122);      //Lambda


    AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);

    AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
    PosAntiv0Daug->SetCutCharge(1);
    AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    NegAntiv0Daug->SetCutCharge(-1);

    Antiv0Cuts->SetCutTransverseRadius(1.4, 200);       // damit v0s erst ab möglichem xiCut berücksichtigt werden - Xi mittlerer flugweg etwa 4,9 cm
    Antiv0Cuts->SetCutCPA(0.97);                        // CPA für sekundäre kleiner?
    Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
    Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
    Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
    Antiv0Cuts->SetPDGCodeNegDaug(2212); //Proton
    Antiv0Cuts->SetPDGCodev0(-3122);     //Lambda

    //  #### Cascade Cuts
    //
    // Xion Cascade
    //
    AliFemtoDreamCascadeCuts *CascadeCutsXion = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
    CascadeCutsXion->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);

    CascadeCutsXion->Setv0Negcuts(XiNegCuts);
    CascadeCutsXion->Setv0PosCuts(XiPosCuts);
    CascadeCutsXion->SetBachCuts(XiBachCuts);
    CascadeCutsXion->SetPDGCodeCasc(3312);      // Xi -
    CascadeCutsXion->SetPDGCodev0(3122);        // Xi 0
    CascadeCutsXion->SetPDGCodePosDaug(2212);   // p +
    CascadeCutsXion->SetPDGCodeNegDaug(-211);   // pi -
    CascadeCutsXion->SetPDGCodeBach(-211);      // pi -

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
    AntiCascadeCutsXion->SetPDGCodeCasc(-3312);     // Xi bar +
    AntiCascadeCutsXion->SetPDGCodev0(-3122);       // Xi bar 0
    AntiCascadeCutsXion->SetPDGCodePosDaug(211);    // pi + 
    AntiCascadeCutsXion->SetPDGCodeNegDaug(-2212);  // p bar -
    AntiCascadeCutsXion->SetPDGCodeBach(211);       // pi +

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
    
    // pairs:   
    // pp                0
    // p bar p           1
    // p Xi              2
    // p bar Xi          3
    // bar p bar p       4
    // bar p Xi          5
    // bar p bar Xi      6
    // Xi Xi             7
    // Xi bar Xi         8
    // bar Xi bar Xi     9
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
    task->Setv0Cuts(v0Cuts);
    task->SetAntiv0Cuts(Antiv0Cuts);
    task->SetTrackCutsXion(CascadeCutsXion);
    task->SetTrackCutsAntiXion(AntiCascadeCutsXion);
    task->SetCollectionConfig(config);

    mgr->AddTask(task);

    TString file = AliAnalysisManager::GetCommonFileName();

    // AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    //
    // mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    //#######
    // ANALYSIS DATA CONTAINERS
    //#######

    // for real particles   -   naming convention for gentle femto : *EvtCuts* - *TrackCuts* - *AntiTrackCuts* - *CascadeCuts* - AntiCascadeCuts* - *Results* - *ResultsQA*
    AliAnalysisDataContainer *coutputEventCuts;
    AliAnalysisDataContainer *coutputProtons;
    AliAnalysisDataContainer *coutputAntiProtons;
    AliAnalysisDataContainer *coutputV0Cuts;
    AliAnalysisDataContainer *coutputAntiV0Cuts;
    AliAnalysisDataContainer *coutputXis;
    AliAnalysisDataContainer *coutputAntiXis;
    AliAnalysisDataContainer *coutputResults;
    AliAnalysisDataContainer *coutputResultsQA;
    coutputEventCuts =      mgr->CreateContainer(Form("EvtCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "EvtCuts"));
    coutputProtons =        mgr->CreateContainer(Form("ProtonTrackCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonTrackCuts"));
    coutputAntiProtons =    mgr->CreateContainer(Form("ProtonAntiTrackCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonAntiTrackCuts"));
    coutputV0Cuts =         mgr->CreateContainer(Form("V0Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("V0Cuts")));
    coutputAntiV0Cuts =     mgr->CreateContainer(Form("AntiV0Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "AntiV0Cuts"));
    coutputXis =            mgr->CreateContainer(Form("XiCascadeCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("XiCascadeCuts")));
    coutputAntiXis =        mgr->CreateContainer(Form("XiAntiCascadeCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("XiAntiCascadeCuts")));
    coutputResults =        mgr->CreateContainer(Form("Results"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "Results"));
    coutputResultsQA =      mgr->CreateContainer(Form("ResultsQA"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ResultsQA"));
    mgr->ConnectOutput(task, 1, coutputEventCuts);
    mgr->ConnectOutput(task, 2, coutputProtons);
    mgr->ConnectOutput(task, 3, coutputAntiProtons);
    mgr->ConnectOutput(task, 4, coutputV0Cuts);
    mgr->ConnectOutput(task, 5, coutputAntiV0Cuts);
    mgr->ConnectOutput(task, 6, coutputXis);
    mgr->ConnectOutput(task, 7, coutputAntiXis);
    mgr->ConnectOutput(task, 8, coutputResults);
    mgr->ConnectOutput(task, 9, coutputResultsQA);

    // for MC   -   naming convention for gentle femto : *TrkCutsMC* - *AntiTrkCutsMC* - *CascCutsMC* - *AntiCascCutsMC*
    if (isMC)
    {
        AliAnalysisDataContainer *coutputTrkCutsMC;
        AliAnalysisDataContainer *coutputAntiTrkCutsMC;
        AliAnalysisDataContainer *coutputCascCutsMC;
        AliAnalysisDataContainer *coutputAntiCascCutsMC;
        coutputTrkCutsMC =      mgr->CreateContainer(Form("ProtonTrkCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonTrkCutsMC"));
        coutputAntiTrkCutsMC =  mgr->CreateContainer(Form("ProtonsAntiTrkCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonsAntiTrkCutsMC"));   
        coutputCascCutsMC =     mgr->CreateContainer(Form("V0CascCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "V0CascCutsMC"));    
        coutputAntiCascCutsMC = mgr->CreateContainer(Form("V0AntiCascCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "V0AntiCascCutsMC"));
        mgr->ConnectOutput(task, 10, coutputTrkCutsMC);
        mgr->ConnectOutput(task, 11, coutputAntiTrkCutsMC);
        mgr->ConnectOutput(task, 12, coutputCascCutsMC);
        mgr->ConnectOutput(task, 13, coutputAntiCascCutsMC);
    }

    return task;
}
