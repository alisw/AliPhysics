#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "AliAnalysisTaskDeuteronCoalescence.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif

//_________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronCoalescence *addTaskDeuteronsCoalescence (Double_t Nch_Transv=7.295)  {
    
    //Get Reweighting
    TFile *inputfile = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/reshaping_protons_pythia/reshaping_protons.root");
    TH1D *hProtWeights = (TH1D*)inputfile->Get("hDataToPythia");
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":deuteron_coalescence";
    
    //Analysis Task
    AliAnalysisTaskDeuteronCoalescence *task = new AliAnalysisTaskDeuteronCoalescence ("task_deuteron_coalescence");
    task -> AliAnalysisTaskDeuteronCoalescence::SetAverageTransverseMultiplicity(Nch_Transv);
    task -> AliAnalysisTaskDeuteronCoalescence::SetReshapingProtons (hProtWeights);
    mgr  -> AddTask(task);
    mgr  -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr  -> ConnectOutput(task,1,mgr->CreateContainer("Coalescence",  TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr  -> ConnectOutput(task,2,mgr->CreateContainer("QAHistograms", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    return task;
}
//_________________________________________________________________________________________________________________________________________________

