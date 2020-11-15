#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "AliAnalysisTaskDeuteronsRT.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif


//____________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronsRT *addTaskDeuteronsRT (
                                                Bool_t    isMC = kFALSE,
                                                UInt_t    triggerType = (AliVEvent::kINT7),//(AliVEvent::kINT7|AliVEvent::kHighMultV0)
                                                Double_t  multMin = 0.0,
                                                Double_t  multMax = 100.0,
                                                Double_t  pt_min_leading = 3.0,
                                                Double_t  average_Nch_Transv = 6.82)  {
    
    
    //Get Analysis Parameters
    TFile *inputFile = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/DEUTERON_ANALYSIS_CUTS/AnalysisParameters_Deuteron_RT.root");
    TH2F *hAnalysisParameters = (TH2F*) inputFile -> Get ("hAnalysisParameters");
    
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
   
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
     
    
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    if ( isMC) fileName += ":deuteron_RT_MC";
    if (!isMC) fileName += ":deuteron_RT_Data";

    
    //Analysis Task
    AliAnalysisTaskDeuteronsRT *task = new AliAnalysisTaskDeuteronsRT ("task_deuterons_vs_RT");
    task -> SelectCollisionCandidates (triggerType);
    task -> AliAnalysisTaskDeuteronsRT::SetTriggerType            (triggerType);
    task -> AliAnalysisTaskDeuteronsRT::SetMultiplicityInterval   (multMin,multMax);
    task -> AliAnalysisTaskDeuteronsRT::SetAverageTransverseMult  (average_Nch_Transv);
    task -> AliAnalysisTaskDeuteronsRT::SetAnalysisParametersSyst (hAnalysisParameters);
    task -> AliAnalysisTaskDeuteronsRT::SetInputData (isMC);
    task -> AliAnalysisTaskDeuteronsRT::SetMinPtLeadingTrack (pt_min_leading);
    mgr  -> AddTask(task);
    mgr  -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr  -> ConnectOutput(task,1,mgr->CreateContainer("Results", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr  -> ConnectOutput(task,2,mgr->CreateContainer("QAplots", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
//____________________________________________________________________________________________________________________________________________

