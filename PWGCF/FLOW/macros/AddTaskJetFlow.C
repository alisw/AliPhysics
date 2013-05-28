///////////////////////////////////////////////////////////////////
//                                                               //            
//                     AddTaskJetFlow                            //
// Author: Redmer A. Bertens, Utrecht University, 2013           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliFlowTrackCuts;

void AddTaskJetFlow( TString name       = "name",
                     TString jets       = "jets",
                     Float_t ptbump     = 0,
                     TArrayI* cent      = 0x0,
                     Float_t MinPOIPt   = 0.15,
                     Float_t MaxPOIPt   = 150,
                     Float_t CCBinsInPt = 100,
                     Bool_t VParticle   = kFALSE,
                     Bool_t debug       = kFALSE  )
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
        Int_t c[] = {0, 10, 30, 50, 70, 90};
        cent = new TArrayI(sizeof(c)/sizeof(c[0]), c);
        printf(" > Setup default centrality binning with %i bins \n ", cent->GetSize());
    }
    // create the cut objects
    AliFlowTrackCuts* CutsRP_VZERO = new AliFlowTrackCuts("CutsRP_VZERO");
    CutsRP_VZERO = CutsRP_VZERO->GetStandardVZEROOnlyTrackCuts();
    AliFlowTrackCuts* CutsRP_TPC = new AliFlowTrackCuts("CutsRP_TPC");
    CutsRP_TPC = CutsRP_TPC->GetStandardTPCStandaloneTrackCuts();
    CutsRP_TPC->SetAODfilterBit(1);
    AliFlowTrackCuts* CutsNull = new AliFlowTrackCuts("CutsNull");
    CutsNull->SetParamType(AliFlowTrackCuts::kGlobal);
    CutsNull->SetEtaRange(+1, -1);
    CutsNull->SetPtRange(+1, -1);
    // add the tasks in a loop, one task for each centrality bin
    for(Int_t i(0); i < cent->GetSize()-1; i++) {
        TString tempName(Form("%s_%i_%i", name.Data(), cent->At(i), cent->At(i+1)));
        // create the task
        AliAnalysisTaskJetFlow* task = new AliAnalysisTaskJetFlow(tempName.Data());   
        task->SetCCMaxPt(MaxPOIPt);
        task->SetCCBinsInPt(CCBinsInPt);
        task->SetDoVParticleAnalysis(VParticle);
        task->SetMinMaxPOIPt(MinPOIPt, MaxPOIPt);
        (debug) ? task->SetDebugMode(1) : task->SetDebugMode(-1);
        if(!task) {
             if(debug) cout << " --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!) " << endl;
             return 0x0;
        }
        else printf(" > added task with name %s and jet collection %s < \n", tempName.Data(), jets.Data());
        task->SetJetCollectionName(jets.Data());
        // pass specific objects and settigns to the task
        task->SetCutsRP(CutsRP_TPC, CutsRP_VZERO);
        task->SetCutsNull(CutsNull);
        task->SetMinMaxCentrality(cent->At(i), cent->At(1+i));
        task->SetPtBump(ptbump);
        mgr->AddTask(task);
        mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
        mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_histograms", tempName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        // connect flow anaysis task
        AliAnalysisDataContainer *flowEvent_VZERO = mgr->CreateContainer(Form("flowEvent_VZERO_%s", tempName.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
        mgr->ConnectOutput(task, 2, flowEvent_VZERO);
        AliAnalysisDataContainer *flowEvent_TPC = mgr->CreateContainer(Form("flowEvent_TPC_%s", tempName.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
        mgr->ConnectOutput(task, 3, flowEvent_TPC);
        TaskJetFlow::AddSPmethod(Form("SPVZERO_A_%s", tempName.Data()), "Qa", 2, flowEvent_VZERO);
        TaskJetFlow::AddSPmethod(Form("SPVZERO_B_%s", tempName.Data()), "Qb", 2, flowEvent_VZERO);
        TaskJetFlow::AddQCmethod(Form("QC_%s", tempName.Data()), 2, flowEvent_TPC);
    }
}
//_____________________________________________________________________________
namespace TaskJetFlow{ // use unique namespace to avoid problems in TRAIN analysis
    void AddSPmethod(char *name, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent)
    {
        TString fileName = AliAnalysisManager::GetCommonFileName();
        fileName+=":SP";
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
        AliAnalysisDataContainer *outSP = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
        AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s", name), kFALSE);
        tskSP->SetApplyCorrectionForNUA(kTRUE);
        tskSP->SetHarmonic(harmonic);
        if(!strcmp("Qa", Qvector)) tskSP->SetTotalQvector("Qa");
        if(!strcmp("Qb", Qvector)) tskSP->SetTotalQvector("Qb");
        tskSP->SetBookOnlyBasicCCH(kFALSE);
        mgr->AddTask(tskSP);
        mgr->ConnectInput(tskSP, 0, flowEvent);
        mgr->ConnectOutput(tskSP, 1, outSP);
    }
//_____________________________________________________________________________
    void AddQCmethod(char *name, int harmonic, AliAnalysisDataContainer *flowEvent)
    {
       TString fileName = AliAnalysisManager::GetCommonFileName();
       fileName+=":QC";
       AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
       AliAnalysisDataContainer *outQC = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
       AliAnalysisTaskQCumulants *tskQC = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s", name), kFALSE);
       tskQC->SetApplyCorrectionForNUA(kTRUE);
       tskQC->SetHarmonic(harmonic);
       tskQC->SetBookOnlyBasicCCH(kFALSE);
       mgr->AddTask(tskQC);
       mgr->ConnectInput(tskQC, 0, flowEvent);
       mgr->ConnectOutput(tskQC, 1, outQC);
    }
//_____________________________________________________________________________
}// end of namespace TaskJetFlow
