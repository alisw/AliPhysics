#ifndef __CINT__
#include <TString.h>
#include <ANALYSIS/AliAnalysisManager.h>
#include <AliRsnAnalysisTask.h>
#include <PWG2/RESONANCES/AliRsnDaughterSelector.h>
#include <PWG2/RESONANCES/AliRsnMiniAnalysisTask.h>
#endif
AliAnalysisTaskSE *AddRsnAnalysisTask(TString format = "esd", Bool_t useMC = kFALSE,Bool_t isMixing = kFALSE,AliRsnInputHandler *rsnIH=0,TList *listRsn=0,TString postfix="")
{
   // create manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error [AddRsnAnalysisTask%s] : mgr is null !!!",postfix.Data()); return 0;}

   // initialize task with all available slots, even if not all of them will be used:
   AliAnalysisTaskSE *task = 0;

   if (gRsnUseMiniPackage) {
      postfix.Prepend("Mini");
      AliRsnMiniAnalysisTask *taskRsnMini = new AliRsnMiniAnalysisTask(Form("Rsn%s",postfix.Data()),useMC);
      taskRsnMini->SetMixPrintRefresh(gRsnMixPrintRefresh);
      Printf("sssssssss %d",gRsnMixPrintRefresh);
      task = (AliAnalysisTaskSE *) taskRsnMini;
   }
   else {
      AliRsnAnalysisTask *taskRsn = new AliRsnAnalysisTask(Form("Rsn%s",postfix.Data()));
      task = (AliAnalysisTaskSE *) taskRsn;
   }

   Printf("pointer is %p",task);

//   gROOT->LoadMacro("RsnConfig.C");
   if (!RsnConfig(task,useMC,isMixing,rsnIH,listRsn)) {
      Printf("Error in RsnConfig.C");
      return 0;
   }

   // add the task to manager
   mgr->AddTask(task);

   AliRsnDaughterSelector *sel = 0;

   if (!gRsnUseMiniPackage) {
      sel = rsnIH->GetSelector();
      sel->Init();
   }

   // connect input container according to source choice
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

   // create paths for the output in the common file
   TString commonPath = AliAnalysisManager::GetCommonFileName();

   // create containers for output
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnHist%s", postfix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath.Data());

   // big outout (expert only)
   Int_t useBigOutput = 0;
   if (useBigOutput) {
      Printf("Using Big output ...");
      task->UseBigOutput();
      output->SetSpecialOutput();
      mgr->SetSpecialOutputLocation("root://lx000.saske.sk:21094//tmp/mvala/");
   }

   mgr->ConnectOutput(task, 1, output);

   return task;
}
