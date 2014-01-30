///////////////////////////////////////////////////////////////////
//                                                               //            
//                     AddTaskJetFlow                            //
// Author: Redmer A. Bertens, Utrecht University, 2013           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliAnalysisTaskRhoVnModulation;

// AddTask macro for the jet flow analysis task
// this task uses an instance AliAnalysisTaskRhoVnModulation
// as a 'master' object for analysis flags and track and event cuts

AliAnalysisTaskRhoVnModulation* AddTaskJetFlow( 
        TString name                       = "name",
        AliAnalysisTaskRhoVnModulation* t  = 0x0,
        Float_t CCMinPt                    = 1,
        Float_t CCMaxPt                    = 150,
        Float_t CCBinsInPt                 = 100,
        Bool_t VParticle                   = kFALSE,
        TArrayD* ptArray                   = 0x0,
        Bool_t VZEROEP                     = kTRUE,
        Bool_t GQC2                        = kTRUE,
        Bool_t QC2                         = kTRUE,
        Bool_t QC4                         = kFALSE,
        Bool_t FlowPackageSP               = kFALSE,
        Bool_t FlowPackageQC               = kTRUE,
        TString LocalRhoName               = "",
        Bool_t debug                       = kFALSE  )
{
    // first check the environment (common to all tasks)
    if(debug) printf("\n\n  >> AddTaskJetFlow <<\n");
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":";
    fileName += t->GetName();
    fileName += "_PWGCF";
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
    AliAnalysisTaskRhoVnModulation* rhoTask = t; 
    // check the centrality setup
    TArrayI* cent(rhoTask->GetCentralityClasses());
    // add the tasks in a loop, one task for each centrality bin
    for(Int_t i(0); i < cent->GetSize()-1; i++) {
        TString tempName(Form("%s_%i_%i", name.Data(), cent->At(i), cent->At(i+1)));
        // create the task
        AliAnalysisTaskJetFlow* task = new AliAnalysisTaskJetFlow(
            tempName.Data(),
            rhoTask,
            VParticle,
            VZEROEP,
            GQC2,            
            QC2,
            QC4,
            FlowPackageSP,
            FlowPackageQC);
        (debug) ? task->SetDebugMode(1) : task->SetDebugMode(-1);
        if(!task) {
             if(debug) printf(" --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!)\n ");
             return 0x0;
        }
        task->SelectCollisionCandidates(rhoTask->GetCollisionCandidates()); 
        task->SetCCMinPt(CCMinPt);
        task->SetCCMaxPt(CCMaxPt);
        task->SetCCBinsInPt(CCBinsInPt);
        task->SetPtBins(ptArray);       // if passed as NULL default a sane default is provided
        task->SetLocalRhoName(LocalRhoName);
        task->SetJetRadius(rhoTask->GetJetRadius());
        // pass specific objects and settigns to the task
        task->SetMinMaxCentrality(cent->At(i), cent->At(1+i));
        mgr->AddTask(task);
        mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
        mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_histograms", tempName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        // connect flow anaysis task
        Bool_t slotTwoFilled(kFALSE);
        if(FlowPackageSP) {
            AliAnalysisDataContainer *flowEvent_VZERO = mgr->CreateContainer(Form("flowEvent_VZERO_%s", tempName.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
            mgr->ConnectOutput(task, 2, flowEvent_VZERO);
            slotTwoFilled = kTRUE;
        }
        if(FlowPackageQC) {
            AliAnalysisDataContainer *flowEvent_TPC = mgr->CreateContainer(Form("flowEvent_TPC_%s", tempName.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
            (slotTwoFilled) ? mgr->ConnectOutput(task, 3, flowEvent_TPC) : mgr->ConnectOutput(task, 2, flowEvent_TPC);
        }
        if(FlowPackageSP) TaskJetFlow::AddSPmethod(Form("SPVZERO_A_%s", tempName.Data()), "Qa", 2, flowEvent_VZERO, t);
        if(FlowPackageSP) TaskJetFlow::AddSPmethod(Form("SPVZERO_B_%s", tempName.Data()), "Qb", 2, flowEvent_VZERO, t);
        if(FlowPackageQC) TaskJetFlow::AddQCmethod(Form("QC_%s", tempName.Data()), 2, flowEvent_TPC, t);
    }
    return rhoTask;
}
//_____________________________________________________________________________
namespace TaskJetFlow{ // use unique namespace to avoid problems in TRAIN analysis
    void AddSPmethod(char *name, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, AliAnalysisTaskRhoVnModulation* t)
    {
        TString fileName = AliAnalysisManager::GetCommonFileName();
        fileName+=":";
        fileName+=t->GetName();
        fileName+="_PWGCF";
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
    void AddQCmethod(char *name, int harmonic, AliAnalysisDataContainer *flowEvent, AliAnalysisTaskRhoVnModulation* t)
    {
       TString fileName = AliAnalysisManager::GetCommonFileName();
       fileName+=":";
       fileName+=t->GetName();
       fileName+="_PWGCF";
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
