#ifndef __CINT__
#include <TString.h>
#include <AliAnalysisManager.h>
#include <AliRsnAnalysisTask.h>
#include <AliRsnDaughterSelector.h>
#include <AliRsnMiniAnalysisTask.h>
#endif
AliAnalysisTaskSE *AddRsnTask(TString rsnPart,TString rsnCut,TString postfix="")
{
   // Analysis Manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error [AddRsnTaskTrain%s] : mgr is null !!!",postfix.Data()); return 0;}

   // initialize task with all available slots, even if not all of them will be used:
   AliAnalysisTaskSE *task = 0;

   Bool_t valid;
   Int_t useMC = AliRsnTrainManager::GetGlobalInt("IsMC",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t physSelBit = AliRsnTrainManager::GetGlobalInt("RsnPhysSelFilterBit",valid);

   if (isRsnMini) {
      postfix.Prepend("Mini");
      AliRsnMiniAnalysisTask *taskRsnMini = new AliRsnMiniAnalysisTask(Form("Rsn%s",postfix.Data()),useMC);
      Int_t refreshPrint = AliRsnTrainManager::GetGlobalInt("RsnMixPrintRefresh",valid);
      if (valid) taskRsnMini->SetMixPrintRefresh(refreshPrint);
      task = (AliAnalysisTaskSE *) taskRsnMini;
   }
   else {
      AliRsnAnalysisTask *taskRsn = new AliRsnAnalysisTask(Form("Rsn%s",postfix.Data()));
      task = (AliAnalysisTaskSE *) taskRsn;
   }

   if (physSelBit>=0) task->SelectCollisionCandidates((AliVEvent::EOfflineTriggerTypes)physSelBit);

   AliRsnInputHandler *rsnIH=0;
   // TODO this is tmp hack
   if (!rsnIH) rsnIH = new AliRsnInputHandler();


   TList *listRsn = new TList();
   listRsn->Add(new TNamed(rsnPart.Data(),rsnCut.Data()));

   gROOT->LoadMacro("RsnConfig.C");
   if (!RsnConfig(task,rsnIH,listRsn)) {
      Printf("Error in RsnConfig.C");
      return 0;
   }

   // add the task to manager
   mgr->AddTask(task);

   AliRsnDaughterSelector *sel = 0;

   if (!isRsnMini) {
      sel = rsnIH->GetSelector();
      sel->Init();
   }

   // connect input container according to source choice
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

   // create paths for the output in the common file
   TString commonPath = AliAnalysisManager::GetCommonFileName();

   // create containers for output
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnHist%s", postfix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath.Data());

//   // big outout (expert only)
//   Int_t useBigOutput = 0;
//   if (useBigOutput) {
//      Printf("Using Big output ...");
//      task->UseBigOutput();
//      output->SetSpecialOutput();
//      mgr->SetSpecialOutputLocation("root://lx000.saske.sk:21094//tmp/mvala/");
//   }

   mgr->ConnectOutput(task, 1, output);

   return task;
}
