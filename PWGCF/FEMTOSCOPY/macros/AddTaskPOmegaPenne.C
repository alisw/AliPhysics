#include "AliAnalysisTaskPOmegaPenne.h"
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskPOmegaPenne *AddTaskPOmegaPenne( bool isMC = false, TString CentEst = "kHM", bool bMixing = false, bool bPairCleanInvMass = true)
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
    // evtCuts->SetMultVsCentPlots(true);

    AliFemtoDreamEventCuts *evtCuts2 = AliFemtoDreamEventCuts::StandardCutsRun2();
    evtCuts2->CleanUpMult(false, false, false, true);
    // evtCuts2->SetMultVsCentPlots(true);

    // Track Cuts
    //
    // AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    // TrackCutsProton->SetFilterBit(128);
    // TrackCutsProton->SetCutCharge(1);

    // AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    // TrackCutsAntiProton->SetFilterBit(128);
    // TrackCutsAntiProton->SetCutCharge(-1);

    // Lambda Cuts
    //
    AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);

    AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);

    AliFemtoDreamv0Cuts *v0Cuts2 = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);

    AliFemtoDreamTrackCuts *Posv0Daug2 = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    AliFemtoDreamTrackCuts *Negv0Daug2 = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);

    // AliFemtoDreamv0Cuts *v0Cuts_rec = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);

    // AliFemtoDreamTrackCuts *Posv0Daug_rec = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    // AliFemtoDreamTrackCuts *Negv0Daug_rec = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);

    // v0Cuts->SetCutInvMass(0.006);                       // same mass range as in XiCuts()
    // v0Cuts->SetCutTransverseRadius(1.4, 200);           // damit v0s erst ab möglichem xiCut berücksichtigt werden - Xi mittlerer flugweg etwa 4,9 cm
    //                                                     // 
    // v0Cuts->SetCutCPA(0.97);                            // CPA für sekundäre weiter!
    // v0Cuts->SetCutDCADaugTov0Vtx(1.4);                  // set a little tighter than dimi (1.5); here it's same to bernie

    v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
    v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
    v0Cuts->SetPDGCodePosDaug(2212); //Proton
    v0Cuts->SetPDGCodeNegDaug(211);  //Pion
    v0Cuts->SetPDGCodev0(3122);      //Lambda

    v0Cuts2->SetPosDaugterTrackCuts(Posv0Daug2);
    v0Cuts2->SetNegDaugterTrackCuts(Negv0Daug2);
    v0Cuts2->SetPDGCodePosDaug(2212); //Proton
    v0Cuts2->SetPDGCodeNegDaug(211);  //Pion
    v0Cuts2->SetPDGCodev0(3122);      //Lambda

    // v0Cuts_rec->SetPosDaugterTrackCuts(Posv0Daug_rec);
    // v0Cuts_rec->SetNegDaugterTrackCuts(Negv0Daug_rec);
    // v0Cuts_rec->SetPDGCodePosDaug(2212); //Proton
    // v0Cuts_rec->SetPDGCodeNegDaug(211);  //Pion
    // v0Cuts_rec->SetPDGCodev0(3122);      //Lambda

    AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);           // V0 cuts - keep Lambda

    AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
    PosAntiv0Daug->SetCutCharge(1);
    AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    NegAntiv0Daug->SetCutCharge(-1);

    AliFemtoDreamv0Cuts *Antiv0Cuts2 = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);          // V0 cuts - keep Xi

    AliFemtoDreamTrackCuts *PosAntiv0Daug2 = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
    PosAntiv0Daug2->SetCutCharge(1);
    AliFemtoDreamTrackCuts *NegAntiv0Daug2 = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    NegAntiv0Daug2->SetCutCharge(-1);

    // Antiv0Cuts->SetCutInvMass(0.006);                   // same v0 mass range as in XiCuts()
    // Antiv0Cuts->SetCutTransverseRadius(1.4, 200);       // damit v0s erst ab möglichem xiCut berücksichtigt werden - Xi mittlerer flugweg etwa 4,9 cm
    // Antiv0Cuts->SetCutCPA(0.97);                        // CPA für sekundäre weiter!
    // Antiv0Cuts->SetCutDCADaugTov0Vtx(1.4);              // set a little tighter than dimi (1.5); here it's same to bernie
    
    Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
    Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
    Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
    Antiv0Cuts->SetPDGCodeNegDaug(2212); //Proton
    Antiv0Cuts->SetPDGCodev0(-3122);     //Lambda

    Antiv0Cuts2->SetPosDaugterTrackCuts(PosAntiv0Daug2);
    Antiv0Cuts2->SetNegDaugterTrackCuts(NegAntiv0Daug2);
    Antiv0Cuts2->SetPDGCodePosDaug(211);  //Pion
    Antiv0Cuts2->SetPDGCodeNegDaug(2212); //Proton
    Antiv0Cuts2->SetPDGCodev0(-3122);     //Lambda

    //  #### Cascade Cuts
    //
    // Xi Cascade
    //
    AliFemtoDreamCascadeCuts *CascadeCutsXion = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
    CascadeCutsXion->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    XiNegCuts->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    XiPosCuts->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    XiBachCuts->SetCheckTPCRefit(false);    //for nanos this is already done while prefiltering (but still mandatory for selection...)

    CascadeCutsXion->Setv0Negcuts(XiNegCuts);
    CascadeCutsXion->Setv0PosCuts(XiPosCuts);
    CascadeCutsXion->SetBachCuts(XiBachCuts);
    CascadeCutsXion->SetPDGCodeCasc(3312);      // Xi -
    CascadeCutsXion->SetPDGCodev0(3122);        // Xi 0
    CascadeCutsXion->SetPDGCodePosDaug(2212);   // p +
    CascadeCutsXion->SetPDGCodeNegDaug(-211);   // pi -
    CascadeCutsXion->SetPDGCodeBach(-211);      // pi -

    // Anti Xi Cascade
    //
    AliFemtoDreamCascadeCuts *AntiCascadeCutsXion = AliFemtoDreamCascadeCuts::XiCuts(isMC, false); 
    AntiCascadeCutsXion->SetXiCharge(1);

    AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    AntiXiNegCuts->SetCutCharge(-1);
    AntiXiNegCuts->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    
    AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    AntiXiPosCuts->SetCutCharge(1);
    AntiXiPosCuts->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    
    AliFemtoDreamTrackCuts *AntiXiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    AntiXiBachCuts->SetCutCharge(1);
    AntiXiBachCuts->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)

    AntiCascadeCutsXion->Setv0Negcuts(AntiXiNegCuts);
    AntiCascadeCutsXion->Setv0PosCuts(AntiXiPosCuts);
    AntiCascadeCutsXion->SetBachCuts(AntiXiBachCuts);
    AntiCascadeCutsXion->SetPDGCodeCasc(-3312);     // Xi bar +
    AntiCascadeCutsXion->SetPDGCodev0(-3122);       // Xi bar 0
    AntiCascadeCutsXion->SetPDGCodePosDaug(211);    // pi + 
    AntiCascadeCutsXion->SetPDGCodeNegDaug(-2212);  // p bar -
    AntiCascadeCutsXion->SetPDGCodeBach(211);       // pi +

    // # 2
    AliFemtoDreamCascadeCuts *CascadeCutsXion2 = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
    CascadeCutsXion2->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts2 = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    XiNegCuts2->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    AliFemtoDreamTrackCuts *XiPosCuts2 = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    XiPosCuts2->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    AliFemtoDreamTrackCuts *XiBachCuts2 = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    XiBachCuts2->SetCheckTPCRefit(false);    //for nanos this is already done while prefiltering (but still mandatory for selection...)

    CascadeCutsXion2->Setv0Negcuts(XiNegCuts2);
    CascadeCutsXion2->Setv0PosCuts(XiPosCuts2);
    CascadeCutsXion2->SetBachCuts(XiBachCuts2);
    CascadeCutsXion2->SetPDGCodeCasc(3312);      // Xi -
    CascadeCutsXion2->SetPDGCodev0(3122);        // Xi 0
    CascadeCutsXion2->SetPDGCodePosDaug(2212);   // p +
    CascadeCutsXion2->SetPDGCodeNegDaug(-211);   // pi -
    CascadeCutsXion2->SetPDGCodeBach(-211);      // pi -

    
    AliFemtoDreamCascadeCuts *AntiCascadeCutsXion2 = AliFemtoDreamCascadeCuts::XiCuts(isMC, false); 
    AntiCascadeCutsXion2->SetXiCharge(1);

    AliFemtoDreamTrackCuts *AntiXiNegCuts2 = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
    AntiXiNegCuts2->SetCutCharge(-1);
    AntiXiNegCuts2->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    
    AliFemtoDreamTrackCuts *AntiXiPosCuts2 = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
    AntiXiPosCuts2->SetCutCharge(1);
    AntiXiPosCuts2->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)
    
    AliFemtoDreamTrackCuts *AntiXiBachCuts2 = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
    AntiXiBachCuts2->SetCutCharge(1);
    AntiXiBachCuts2->SetCheckTPCRefit(false);     //for nanos this is already done while prefiltering (but still mandatory for selection...)

    AntiCascadeCutsXion2->Setv0Negcuts(AntiXiNegCuts2);
    AntiCascadeCutsXion2->Setv0PosCuts(AntiXiPosCuts2);
    AntiCascadeCutsXion2->SetBachCuts(AntiXiBachCuts2);
    AntiCascadeCutsXion2->SetPDGCodeCasc(-3312);     // Xi bar +
    AntiCascadeCutsXion2->SetPDGCodev0(-3122);       // Xi bar 0
    AntiCascadeCutsXion2->SetPDGCodePosDaug(211);    // pi + 
    AntiCascadeCutsXion2->SetPDGCodeNegDaug(-2212);  // p bar -
    AntiCascadeCutsXion2->SetPDGCodeBach(211);       // pi +

    std::vector<int> PDGParticles;
    // PDGParticles.push_back(2212);   // Protons
    // PDGParticles.push_back(2212);

    PDGParticles.push_back(3122);   // Lamdas
    PDGParticles.push_back(3122);

    PDGParticles.push_back(3312);   // Xis
    PDGParticles.push_back(3312);

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

    int iQaPairs = 10;
    for (int i = 0; i < iQaPairs; ++i)
    {
        pairQA.push_back(0);
        closeRejection.push_back(false);
        NBins.push_back(750);
        kMin.push_back(0.);
        kMax.push_back(3.);
    }
    // // // // pairs:   
    // // // // pp                0
    // // // // p bar p           1
    // // // // p Xi              2
    // // // // p bar Xi          3
    // // // // bar p bar p       4
    // // // // bar p Xi          5
    // // // // bar p bar Xi      6
    // // // // Xi Xi             7
    // // // // Xi bar Xi         8
    // // // // bar Xi bar Xi     9
    
    // pairQA[0] = 11;               // protons (old)
    // pairQA[4] = 11; 
    // pairQA[2] = 13;     
    // pairQA[6] = 13;    
    // closeRejection[0] = true;
    // closeRejection[4] = true;     // ##

    // pairs:   

    // LL                0
    // L bar L           1
    // L Xi              2
    // L bar Xi          3
    // bar L bar L       4
    // bar L Xi          5
    // bar L bar Xi      6
    // Xi Xi             7
    // Xi bar Xi         8
    // bar Xi bar Xi     9

    pairQA[0] = 22;     // LL                0
    pairQA[1] = 22;     // L bar L           1
    pairQA[4] = 22;     // bar L bar L       4

    pairQA[2] = 23;     // L Xi              2
    pairQA[3] = 23;     // L bar Xi          3
    pairQA[5] = 23;     // bar L Xi          5
    pairQA[6] = 23;     // bar L bar Xi      6
    
    pairQA[7] = 33;     // Xi Xi             7
    pairQA[8] = 33;     // Xi bar Xi         8
    pairQA[9] = 33;     // bar Xi bar Xi     9
        
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

    // full blast QA
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMassQA(false);
    // config->SetMinimalBookingME(true);

    // ##### Task creation!!!!! ################
    AliAnalysisTaskPOmegaPenne *task = new AliAnalysisTaskPOmegaPenne("FemtoDreamPOmegaPenne", isMC, bMixing, bPairCleanInvMass);
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
    task->SetTrackCutsXion(CascadeCutsXion);
    task->SetTrackCutsAntiXion(AntiCascadeCutsXion);
    // # 2
    task->SetEventCuts2(evtCuts2);
    task->Setv0Cuts2(v0Cuts2);
    task->SetAntiv0Cuts2(Antiv0Cuts2);
    task->SetTrackCutsXion2(CascadeCutsXion2);
    task->SetTrackCutsAntiXion2(AntiCascadeCutsXion2);
    // # recombination cuts
    // task->Setv0Cuts_rec(v0Cuts_rec);
    // task->SetAntiv0Cuts_rec(Antiv0Cuts_rec);

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
    // for real particles   -   naming convention for gentle femto : *EvtCuts* - *TrackCuts* - *AntiTrackCuts* - *CascadeCuts* - *AntiCascadeCuts* - *Results* - *ResultsQA*
    AliAnalysisDataContainer *coutputEventCuts;
    // AliAnalysisDataContainer *coutputProtons;
    // AliAnalysisDataContainer *coutputAntiProtons;
    AliAnalysisDataContainer *coutputV0Cuts;
    AliAnalysisDataContainer *coutputAntiV0Cuts;
    AliAnalysisDataContainer *coutputXis;
    AliAnalysisDataContainer *coutputAntiXis;
    AliAnalysisDataContainer *coutputResults;
    AliAnalysisDataContainer *coutputResultsQA;
    // #2 only keep Lambda instead of Xi
    AliAnalysisDataContainer *coutputEventCuts2;
    AliAnalysisDataContainer *coutputV0Cuts2;
    AliAnalysisDataContainer *coutputAntiV0Cuts2;
    AliAnalysisDataContainer *coutputXis2;
    AliAnalysisDataContainer *coutputAntiXis2;
    AliAnalysisDataContainer *coutputResultsQA2;
    AliAnalysisDataContainer *coutputResults2;
    
    AliAnalysisDataContainer *coutputRecombBeforePairclean;     // recombination statistics BEFORE PairCleaner
    AliAnalysisDataContainer *coutputRecombAfterPairclean;       // recombination statistics AFTER PairCleaner

    coutputEventCuts =      mgr->CreateContainer(Form("EvtCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "EvtCuts"));
    // coutputProtons =        mgr->CreateContainer(Form("ProtonTrackCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonTrackCuts"));
    // coutputAntiProtons =    mgr->CreateContainer(Form("ProtonAntiTrackCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonAntiTrackCuts"));
    coutputV0Cuts =         mgr->CreateContainer(Form("V0Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("V0Cuts")));
    coutputAntiV0Cuts =     mgr->CreateContainer(Form("AntiV0Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "AntiV0Cuts"));
    coutputXis =            mgr->CreateContainer(Form("XiCascadeCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("XiCascadeCuts")));
    coutputAntiXis =        mgr->CreateContainer(Form("XiAntiCascadeCuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("XiAntiCascadeCuts")));
    coutputResults =        mgr->CreateContainer(Form("Results"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "Results"));
    coutputResultsQA =      mgr->CreateContainer(Form("ResultsQA"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ResultsQA"));
    
    coutputEventCuts2 =      mgr->CreateContainer(Form("EvtCuts2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "EvtCuts2"));
    coutputV0Cuts2 =         mgr->CreateContainer(Form("V0Cuts2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("V0Cuts2")));
    coutputAntiV0Cuts2 =     mgr->CreateContainer(Form("AntiV0Cuts2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "AntiV0Cuts"));
    coutputXis2 =            mgr->CreateContainer(Form("XiCascadeCuts2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("XiCascadeCuts2")));
    coutputAntiXis2 =        mgr->CreateContainer(Form("XiAntiCascadeCuts2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), Form("XiAntiCascadeCuts2")));
    coutputResultsQA2 =      mgr->CreateContainer(Form("ResultsQA2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ResultsQA2"));
    coutputResults2 =        mgr->CreateContainer(Form("Results2"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "Results2"));
    
    coutputRecombBeforePairclean =   mgr->CreateContainer(Form("RecombinationBeforePairClean"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "RecombinationBeforePairClean"));
    coutputRecombAfterPairclean =   mgr->CreateContainer(Form("RecombinationAfterPairClean"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "RecombinationAfterPairClean"));

    mgr->ConnectOutput(task, 1, coutputEventCuts);
    // mgr->ConnectOutput(task, 2, coutputProtons);
    // mgr->ConnectOutput(task, 3, coutputAntiProtons);
    mgr->ConnectOutput(task, 2, coutputV0Cuts);
    mgr->ConnectOutput(task, 3, coutputAntiV0Cuts);
    mgr->ConnectOutput(task, 4, coutputXis);
    mgr->ConnectOutput(task, 5, coutputAntiXis);
    mgr->ConnectOutput(task, 6, coutputResults);
    mgr->ConnectOutput(task, 7, coutputResultsQA);    // paircleaner - keep lambda not Xi
    mgr->ConnectOutput(task, 8, coutputEventCuts2);
    mgr->ConnectOutput(task, 9, coutputV0Cuts2);
    mgr->ConnectOutput(task, 10, coutputAntiV0Cuts2);
    mgr->ConnectOutput(task, 11, coutputXis2);
    mgr->ConnectOutput(task, 12, coutputAntiXis2);
    mgr->ConnectOutput(task, 13, coutputResults2);
    mgr->ConnectOutput(task, 14, coutputResultsQA2);    // paircleaner - keep Xi not lambda
    mgr->ConnectOutput(task, 15, coutputRecombBeforePairclean);    // recombination statistics
    mgr->ConnectOutput(task, 16, coutputRecombAfterPairclean);    // recombination statistics
    

    // for MC   -   naming convention for gentle femto : *TrkCutsMC* - *AntiTrkCutsMC* - *CascCutsMC* - *AntiCascCutsMC*
    if (isMC)
    {
        // AliAnalysisDataContainer *coutputTrkCutsMC;
        // AliAnalysisDataContainer *coutputAntiTrkCutsMC;
        AliAnalysisDataContainer *coutputCascCutsMC;
        AliAnalysisDataContainer *coutputAntiCascCutsMC;
        // coutputTrkCutsMC =      mgr->CreateContainer(Form("ProtonTrkCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonTrkCutsMC"));
        // coutputAntiTrkCutsMC =  mgr->CreateContainer(Form("ProtonsAntiTrkCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "ProtonsAntiTrkCutsMC"));   
        coutputCascCutsMC =     mgr->CreateContainer(Form("V0CascCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "V0CascCutsMC"));    
        coutputAntiCascCutsMC = mgr->CreateContainer(Form("V0AntiCascCutsMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), "V0AntiCascCutsMC"));
        // mgr->ConnectOutput(task, 10, coutputTrkCutsMC);
        // mgr->ConnectOutput(task, 11, coutputAntiTrkCutsMC);
        mgr->ConnectOutput(task, 17, coutputCascCutsMC);
        mgr->ConnectOutput(task, 18, coutputAntiCascCutsMC);
    }
    return task;
}
