#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "AliAnalysisTaskXiNucleusInteraction.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif

//_________________________________________________________________________________________________________________________________________________
AliAnalysisTaskXiNucleusInteraction *addTaskXiNucleusInteraction (UInt_t trigger=(AliVEvent::kINT7),
                                                                  Double_t multLow=0.0,
                                                                  Double_t multHigh=100.0)  {
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":YN_Interaction";
    
    //Analysis Task
    AliAnalysisTaskXiNucleusInteraction *task = new AliAnalysisTaskXiNucleusInteraction ("xi_nucleus_interaction_task");
    task -> SelectCollisionCandidates (trigger);
    task -> AliAnalysisTaskXiNucleusInteraction::SetMultiplicity(trigger,multLow,multHigh);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("Results", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("QAPlots", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
//_________________________________________________________________________________________________________________________________________________

