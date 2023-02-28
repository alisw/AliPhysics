#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "AliAnalysisTaskSimpleCoalescenceHelium3.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif

//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceHelium3 *addTaskHelium3SimpleCoalescence ()  {
    
    //Get Input File
    TFile *input = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/reshaping_protons_pythia/reshaping_protons.root");
    
    //Get Weights
    TF1 *fProtWeights = (TF1*) input -> Get("protonWeight");
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString filename = AliAnalysisManager::GetCommonFileName();
    filename += ":helium3_coalescence";

    //Analysis Task
    AliAnalysisTaskSimpleCoalescenceHelium3 *task = new AliAnalysisTaskSimpleCoalescenceHelium3 ("task_helium3_coalescence");
    task -> AliAnalysisTaskSimpleCoalescenceHelium3::SetProtonWeights (fProtWeights);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("Results",TList::Class(),AliAnalysisManager::kOutputContainer,filename.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("QAPlots",TList::Class(),AliAnalysisManager::kOutputContainer,filename.Data()));
    
    return task;
}
//_______________________________________________________________________________________________________________________________________

