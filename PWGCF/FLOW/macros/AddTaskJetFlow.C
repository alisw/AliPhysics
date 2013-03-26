///////////////////////////////////////////////////////////////////
//                                                               //            
//                     AddTaskJetFlow                            //
// Author: Redmer A. Bertens, Utrecht University, 2013           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliFlowTrackCuts;

void AddTaskJetFlow( TString name    = "name",
                     TString jets    = "jets",
                     Float_t ptbump  = 200,
                     TArrayI* cent   = 0x0,
                     Bool_t debug    = kTRUE  )
{
    // first check the environment (common to all tasks)
    if(debug) printf("\n\n  >> AddTaskJetFlow <<\n");
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":JetFlow";
    if(debug) printf("      - filename: %s \n",fileName.Data());
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        if(debug) cout << " Fatal error: no analysis manager found! " << endl;
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        if(debug) cout << " Fatal error: no imput event handler found!" << endl;
        return 0x0;
    }
    // check the centrality setup
    if(!cent) {
        Int_t c[] = {0, 20, 40, 60, 80, 100};
        cent = new TArrayI(sizeof(c)/sizeof(c[0]), c);
        printf(" > Setup default centrality binning with %i bins \n ", cent->GetSize());
    }
    // create the cut objects
    AliFlowTrackCuts* CutsRP = new AliFlowTrackCuts("CutsRP");
    CutsRP = CutsRP->GetStandardVZEROOnlyTrackCuts();
    AliFlowTrackCuts* CutsNull = new AliFlowTrackCuts("CutsNull");
    CutsNull->SetParamType(AliFlowTrackCuts::kGlobal);
    CutsNull->SetEtaRange(+1, -1);
    CutsNull->SetPtRange(+1, -1);
    // add the tasks in a loop, one task for each centrality bin
    for(Int_t i(0); i < cent->GetSize()-1; i++) {
        TString tempName(Form("%s_%i_%i", name.Data(), cent->At(i), cent->At(i+1)));
        // create the task
        AliAnalysisTaskJetFlow* task = new AliAnalysisTaskJetFlow(tempName.Data());   
        if(!task) {
             if(debug) cout << " --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!) " << endl;
             return 0x0;
        }
        else printf(" > added task with name %s and jet collection %s < \n", tempName.Data(), jets.Data());
        task->SetJetCollectionName(jets.Data());
        // pass specific objects and settigns to the task
        task->SetCutsRP(CutsRP);
        task->SetCutsNull(CutsNull);
        task->SetMinMaxCentrality(cent->At(i), cent->At(1+i));
        task->SetPtBump(ptbump);
        mgr->AddTask(task);
        mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
        mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_histograms", tempName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     
   // connect flow anaysis task
        AliAnalysisDataContainer *flowEvent = mgr->CreateContainer(Form("flowEvent_%s", tempName.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
        mgr->ConnectOutput(task, 2, flowEvent);
        TaskJetFlow::AddSPmethod(Form("SP_A_%s", tempName.Data()), -0.8, -0, 0, +0.8, "Qa", 2, flowEvent);
        TaskJetFlow::AddSPmethod(Form("SP_B_%s", tempName.Data()), -0.8, -0, 0, +0.8, "Qb", 2, flowEvent);
    }
}
//_____________________________________________________________________________
namespace TaskJetFlow{ // use unique namespace to avoid problems in TRAIN analysis
    void AddSPmethod(char *name, double minEtaA, double maxEtaA, double minEtaB, double maxEtaB, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent)
    {
       // add sp task and filter tasks
       AliFlowTrackSimpleCuts* a = new AliFlowTrackSimpleCuts("a");
       a->SetEtaMin(minEtaA);
       a->SetEtaMax(maxEtaA);
       a->SetMassMin(-999);
       a->SetMassMax(999);
       AliFlowTrackSimpleCuts* b = new AliFlowTrackSimpleCuts("b");
       b->SetEtaMin(minEtaB);
       b->SetEtaMax(maxEtaB);
       b->SetMassMin(-999);
       b->SetMassMax(999);
       TString fileName = AliAnalysisManager::GetCommonFileName();
       fileName+=":SP";
       fileName+="name";
       AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
       AliAnalysisDataContainer *filteredFlowEvent = mgr->CreateContainer(Form("Filter_%s", name), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
       AliAnalysisTaskFilterFE *tskFilter(0x0);
       if(!strcmp("Qa", Qvector)) tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", name), 0x0, a);
       if(!strcmp("Qb", Qvector)) tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", name), 0x0, b);
       if(tskFilter) tskFilter->SetSubeventEtaRange(-10, 0, 0, 10);
       else return;
       mgr->AddTask(tskFilter);
       printf( " > Added FILTER TASK < \n ");
       mgr->ConnectInput(tskFilter, 0, flowEvent);
       mgr->ConnectOutput(tskFilter, 1, filteredFlowEvent);
       AliAnalysisDataContainer *outSP = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
       AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s", name), kFALSE);
       tskSP->SetApplyCorrectionForNUA(kTRUE);
       tskSP->SetHarmonic(harmonic);
       tskSP->SetTotalQvector(Qvector);
       tskSP->SetBookOnlyBasicCCH(kFALSE);
       mgr->AddTask(tskSP);
       mgr->ConnectInput(tskSP, 0, filteredFlowEvent);
       mgr->ConnectOutput(tskSP, 1, outSP);
    }
}// end of namespace TaskJetFlow
