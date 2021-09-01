#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSESemileptonicOmegac0KFP.h"
#include "AliNormalizationCounter.h"
#include <TTree.h>
#include <TString.h>
#include <TList.h>
#include "TFile.h"
#endif

AliAnalysisTaskSESemileptonicOmegac0KFP *AddTaskOmegac2eleOmegafromKFP(TString finname="", Bool_t theMCon=kFALSE, Bool_t writeQATree=kFALSE, TString cuttype="")

{

    Bool_t writeOmegac0RecTree = kTRUE;
    Bool_t writeOmegac0MCGenTree = kFALSE;
    if(theMCon) writeOmegac0MCGenTree = kTRUE;
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      ::Error("AddMyTask", "No analysis manager to connect to.");
      return NULL;
    }
    
   if (!mgr->GetInputEventHandler()) {
      return NULL;
     }

   Bool_t stdcuts=kFALSE;
   TFile* filecuts;
   if (finname.EqualTo("") ) {
     stdcuts=kTRUE;
   } else {
     filecuts=TFile::Open(finname.Data());
     if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
       printf("Input file not found : check your cut object");
     }
   }
    
   AliRDHFCutsOmegactoeleOmegafromKFP * RDHFCutsOmegac2eleOmegaanal = new AliRDHFCutsOmegactoeleOmegafromKFP();
   if (stdcuts) RDHFCutsOmegac2eleOmegaanal->SetStandardCutsPP2010();
   else RDHFCutsOmegac2eleOmegaanal = (AliRDHFCutsOmegactoeleOmegafromKFP*)filecuts->Get("eleOmegaAnalysisCuts");
   RDHFCutsOmegac2eleOmegaanal->SetName("eleOmegaAnalysisCuts");
   RDHFCutsOmegac2eleOmegaanal->SetMinPtCandidate(-1.);
   RDHFCutsOmegac2eleOmegaanal->SetMaxPtCandidate(10000.);
   if (!RDHFCutsOmegac2eleOmegaanal) {
     cout << "Specific AliRDHFCutsOmegac2eleOmegaanal not found\n";
     return NULL;
   }
     printf("^^^^^^^^^^^^CREATE TASK\n");
  
    // create an instance of your task
   AliAnalysisTaskSESemileptonicOmegac0KFP *task = new AliAnalysisTaskSESemileptonicOmegac0KFP("AliAnalysisTaskSESemileptonicOmegac0KFP",RDHFCutsOmegac2eleOmegaanal);
   if(!task) return NULL;
    
   task->SetMC(theMCon);
   task->SetDebugLevel(1);
   task->SetWriteOmegac0MCGenTree(writeOmegac0MCGenTree);
   task->SetWriteOmegac0Tree(writeOmegac0RecTree);
   task->SetWriteOmegac0QATree(writeQATree);
   mgr->AddTask(task);
    

 // Create and connect containers for input/output
   TString outputfile = AliAnalysisManager::GetCommonFileName();
   outputfile += ":PWGHF_D2H_Omegac2eleOmegaKFP_";
   outputfile += cuttype.Data();
       
   mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("CutObj_%s", cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
   mgr->ConnectOutput(task, 2, mgr->CreateContainer(Form("Counter_%s", cuttype.Data()), AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
   mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("AnaHisto_%s", cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
   mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("tree_event_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
   mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("tree_Omegac0_%s",cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
   mgr->ConnectOutput(task,6,mgr->CreateContainer(Form("tree_Omegac0_MCGen_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
   mgr->ConnectOutput(task,7,mgr->CreateContainer(Form("tree_Omegac0_QA_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

 // in the end, this macro returns a pointer to your task. this will be convenient later on
 // when you will run your analysis in an analysis train on grid
   return task;

}
