///////////////////////////////////////////////////////////////////
//                                                               //            
//                     AddTaskJetFlow                            //
// Author: Redmer A. Bertens, Utrecht University, 2013           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliFlowTrackCuts;

AliAnalysisTaskJetFlow* AddTaskJetFlow( TString name = "name",
                                        TString jets = "jets",
                                        Bool_t debug = kTRUE  )

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
 
   // create the task
   AliAnalysisTaskJetFlow* task = new AliAnalysisTaskJetFlow(name.Data());   
   if(!task) {
        if(debug) cout << " --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!) " << endl;
        return 0x0;
    }
   else printf(" > added task with name %s and jet collection %s < \n", name.Data(), jets);
   task->SetJetCollectionName(jets.Data());

   // pass specific objects to the task
   AliFlowTrackCuts* CutsRP = new AliFlowTrackCuts("CutsRP");
   CutsRP = CutsRP->GetStandardVZEROOnlyTrackCuts();
   task->SetCutsRP(CutsRP);

   AliFlowTrackCuts* CutsNull = new AliFlowTrackCuts("CutsNull");
   CutsNull->SetParamType(AliFlowTrackCuts::kGlobal);
   CutsNull->SetEtaRange(+1, -1);
   CutsNull->SetPtRange(+1, -1);
   task->SetCutsNull(CutsNull);
   
   mgr->AddTask(task);
   mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task,1,mgr->CreateContainer("GenericOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

   // connect flow anaysis task
   AliAnalysisDataContainer *flowEvent = mgr->CreateContainer("flowEvent", AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
   mgr->ConnectOutput(task, 2, flowEvent);
   TaskJetFlow::AddSPmethod("SP_A", -0.8, -0, 0, +0.8, "Qa", 2, flowEvent);
   TaskJetFlow::AddSPmethod("SP_B", -0.8, -0, 0, +0.8, "Qb", 2, flowEvent);

   return task;

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
