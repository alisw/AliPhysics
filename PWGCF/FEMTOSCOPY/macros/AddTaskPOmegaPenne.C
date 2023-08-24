#include "AliAnalysisTaskPOmegaPenne.h"
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskPOmegaPenne *AddTaskPOmegaPenne(bool isMC = true, 
                                               TString CentEst = "kHM",
                                               bool mixBeforePC = false,
                                               bool mixAfterPC = true, 
                                               bool isInvMassPairClean = true,
                                               TString  multTrigger = "1",
                                               bool fullBlastQA = true
                                               )
{
    // #define RUN_SECOND_SET_OF_CUTS
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
    // evtCuts->SetMultVsCentPlots(true);

    // Lambda Cuts
    //
    AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);
    // v0Cuts->SetCutInvMass(0.050);
    // v0Cuts->SetAxisInvMassPlots(600, 0.95, 1.25);

    AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);

    v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
    v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
    v0Cuts->SetPDGCodePosDaug(2212); //Proton
    v0Cuts->SetPDGCodeNegDaug(211);  //Pion
    v0Cuts->SetPDGCodev0(3122);      //Lambda


    AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);           // V0 cuts - keep Lambda
    // Antiv0Cuts->SetCutInvMass(0.050);
    // Antiv0Cuts->SetAxisInvMassPlots(600, 0.95, 1.25);

    AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
    PosAntiv0Daug->SetCutCharge(1);
    AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    NegAntiv0Daug->SetCutCharge(-1);

    
    // Antiv0Cuts->SetCutInvMass(0.006);                   // same v0 mass range as in XiCuts()
    // Antiv0Cuts->SetCutTransverseRadius(1.4, 200);       // damit v0s erst ab möglichem xiCut berücksichtigt werden - Xi mittlerer flugweg etwa 4,9 cm
    // Antiv0Cuts->SetCutCPA(0.97);                        // CPA für sekundäre weiter!
    // Antiv0Cuts->SetCutDCADaugTov0Vtx(1.4);              // set a little tighter than dimi (1.5); here it's same to bernie

    Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
    Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
    Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
    Antiv0Cuts->SetPDGCodeNegDaug(2212); //Proton
    Antiv0Cuts->SetPDGCodev0(-3122);     //Lambda

    //////////////////////////////////////
    // Config
    /////////////////////////////////////
    std::vector<int> PDGParticles;

    PDGParticles.push_back(3122);   // Lamdas
    PDGParticles.push_back(3122);

    /* std::vector( size_type count, const T& value, const Allocator& alloc = Allocator()); */
    // std::vector<int> NBins = std::vector<int>(10, 750);
    // std::vector<float> kMin = std::vector<float>(10, 0.);
    // std::vector<float> kMax = std::vector<float>(10, 3.);
    // std::vector<int> pairQA = std::vector<int>(10, 0);
    // std::vector<bool> closeRejection = std::vector<bool>(10, false);

    std::vector<int> NBins;
    std::vector<float> kMin;
    std::vector<float> kMax;
    std::vector<int> pairQA;
    std::vector<bool> closeRejection;

    int iQaPairs = 3;
    for (int i = 0; i < iQaPairs; ++i)
    {
        pairQA.push_back(0);
        closeRejection.push_back(false);
        NBins.push_back(2000);
        kMin.push_back(0.);
        kMax.push_back(2.);
    }
    
    // pairs:   

    // LL                0
    // L bar L           1
    // bar L bar L       2

    pairQA[0] = 22;     // LL                0
    pairQA[1] = 22;     // L bar L           1
    pairQA[2] = 22;     // bar L bar L       2

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
    config->SetUseEventMixing(true);        // defaults to true
    config->SetExtendedQAPairs(pairQA);
    config->SetClosePairRejection(closeRejection);
    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    if (multTrigger == "0")
    {
        config->SetMultiplicityEstimator(AliFemtoDreamEvent::kSPD);
    }
    if (multTrigger == "1")   // default
    {
        config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
    }
    if (multTrigger == "2")
    {
        config->SetMultiplicityEstimator(AliFemtoDreamEvent::kV0M);
    }
    if (multTrigger == "3")
    {
        config->SetMultiplicityEstimator(AliFemtoDreamEvent::kV0A);
    }
    if (multTrigger == "4")
    {
        config->SetMultiplicityEstimator(AliFemtoDreamEvent::kV0C);
    }
    

    if (isMC) 
    {
        config->SetMomentumResolution(true);
    }

    // full blast QA
    if (!fullBlastQA)
    {
        evtCuts->SetMinimalBooking(true);
        v0Cuts->SetMinimalBooking(true);
        Antiv0Cuts->SetMinimalBooking(true);

        config->SetMinimalBookingME(true);
        config->SetMinimalBookingSample(true);
    }
    
    if(fullBlastQA)
    {
        config->SetkTBinning(true);
        config->SetPtQA(true);
        config->SetMassQA(true);
    }

    // ##### Task creation!!!!! ################
    AliAnalysisTaskPOmegaPenne *task = new AliAnalysisTaskPOmegaPenne("FemtoDreamPOmegaPenne", isMC);

    task->SetMixBeforePC(mixBeforePC);
    task->SetMixAfterPC(mixAfterPC);
    task->SetFullBlastQA(fullBlastQA);
    task->SetInvMassPairClean(isInvMassPairClean);
    task->SetMultTrigger(multTrigger);

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
    //#
    task->SetEventCuts(evtCuts);
    // task->SetTrackCutsProton(TrackCutsProton);
    // task->SetTrackCutsAntiProton(TrackCutsAntiProton);
    task->Setv0Cuts(v0Cuts);
    task->SetAntiv0Cuts(Antiv0Cuts);

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
    AliAnalysisDataContainer *coutputV0Cuts;
    AliAnalysisDataContainer *coutputAntiV0Cuts;
    AliAnalysisDataContainer *coutputResults;
    AliAnalysisDataContainer *coutputResultsQA;
    AliAnalysisDataContainer *coutputRecombAfterPairclean;       // recombination statistics AFTER PairCleaner
    AliAnalysisDataContainer *coutputResults2;              
    AliAnalysisDataContainer *coutputResults3;            

// for MC   -   naming convention for gentle femto : *TrkCutsMC* - *AntiTrkCutsMC* - *CascCutsMC* - *AntiCascCutsMC*
    AliAnalysisDataContainer *coutputv0CutsMC;
    AliAnalysisDataContainer *coutputAntiv0CutsMC;

    coutputEventCuts =      mgr->CreateContainer(Form("EvtCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "EvtCuts"));
    coutputV0Cuts =         mgr->CreateContainer(Form("V0Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("V0Cuts")));
    coutputAntiV0Cuts =     mgr->CreateContainer(Form("AntiV0Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "AntiV0Cuts"));
    coutputResults =        mgr->CreateContainer(Form("Results"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "Results"));
    coutputResultsQA =      mgr->CreateContainer(Form("ResultsQA_CPAClean"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ResultsQA_CPAClean"));
    
    coutputRecombAfterPairclean =   mgr->CreateContainer(Form("RecombinationAfterPairClean"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "RecombinationAfterPairClean"));
    
    coutputResults2 =        mgr->CreateContainer(Form("ResultsCleanInvMass"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ResultsCleanInvMass"));
    coutputResults3 =        mgr->CreateContainer(Form("ResultsCleanAtRandom"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ResultsCleanAtRandom"));

    mgr->ConnectOutput(task, 1, coutputEventCuts);
    mgr->ConnectOutput(task, 2, coutputV0Cuts);
    mgr->ConnectOutput(task, 3, coutputAntiV0Cuts);
    mgr->ConnectOutput(task, 4, coutputResults);
    mgr->ConnectOutput(task, 5, coutputResultsQA);
    mgr->ConnectOutput(task, 6, coutputRecombAfterPairclean);
    mgr->ConnectOutput(task, 7, coutputResults2);
    mgr->ConnectOutput(task, 8, coutputResults3);   // recombination statistics
    if (isMC)
    {
        coutputv0CutsMC =     mgr->CreateContainer(Form("v0CutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "v0CutsMC"));    
        coutputAntiv0CutsMC = mgr->CreateContainer(Form("Antiv0CutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "Antiv0CutsMC"));
        mgr->ConnectOutput(task, 9, coutputv0CutsMC);
        mgr->ConnectOutput(task, 10, coutputAntiv0CutsMC);
    }
    
    
    return task;
}
