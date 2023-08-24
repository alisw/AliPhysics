#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "AliAnalysisTaskSimpleCoalescenceDeuteronInJets.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif

//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceDeuteronInJets *addTaskDeuteronsInJets (Double_t jetRadius=0.5, Double_t ptTrigger=5.0)  {
    
    //Get Input File
    TFile *input = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/reshaping_protons_pythia/weights_protons_in_jet.root");
    
    //Get Weights
    TF1 *fProtWeights = (TF1*) input -> Get("fWeight");
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString filename = AliAnalysisManager::GetCommonFileName();
    filename += Form(":jetRadius%.1f_ptMax%.1f",jetRadius,ptTrigger);

    //Analysis Task
    AliAnalysisTaskSimpleCoalescenceDeuteronInJets *task = new AliAnalysisTaskSimpleCoalescenceDeuteronInJets (Form("taskDeutInJets_jetRadius%.1f_ptMax%.1f",jetRadius,ptTrigger));
    task -> AliAnalysisTaskSimpleCoalescenceDeuteronInJets::SetProtonWeights (fProtWeights);
    task -> AliAnalysisTaskSimpleCoalescenceDeuteronInJets::SetJetRadius (jetRadius);
    task -> AliAnalysisTaskSimpleCoalescenceDeuteronInJets::SetMaximumPt (ptTrigger);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer(Form("Results_jetRadius%.1f_ptMax%.1f",jetRadius,ptTrigger),TList::Class(),AliAnalysisManager::kOutputContainer,filename.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer(Form("QAPlots_jetRadius%.1f_ptMax%.1f",jetRadius,ptTrigger),TList::Class(),AliAnalysisManager::kOutputContainer,filename.Data()));
    
    return task;
}
//_______________________________________________________________________________________________________________________________________

