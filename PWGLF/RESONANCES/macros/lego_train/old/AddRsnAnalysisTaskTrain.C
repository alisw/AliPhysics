#ifndef __CINT__
#include <TString.h>
#endif
AliAnalysisTaskSE *AddRsnAnalysisTaskTrain(TString addRsnManager="AddRsnToManagerTrain.C",
                                           Int_t isESD=0,
                                           Int_t isMC=0,
                                           Int_t isPP=1,
                                           Int_t useRsnMini = 1,
                                           Int_t useMixing = 0,
                                           Int_t numMix = 10,
                                           Int_t fullOutput = 1,
                                           AliRsnInputHandler *rsnIH=0,
                                           TString postfix="")
{
   // create manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error [AddRsnAnalysisTask%s] : mgr is null !!!",postfix.Data()); return 0;}

   // initialize task with all available slots, even if not all of them will be used:
   AliAnalysisTaskSE *task = 0;

   TString lego_path = AliAnalysisManager::GetGlobalStr("rsnLegoTrainPath",valid);
   TString rsnManagerMacro = "RsnManagerTrain.C";
   gROOT->LoadMacro(gSystem->EpandPathName(TString::Format("%s/%s",lego_path.Data(),rsnManagerMacro.Data()).Data()));
   TList *listRsn = RsnManagerTrain(addRsnManager,isESD,isMC,isPP,useRsnMini,useMixing,numMix,fullOutput);
   if (!listRsn) return 0;

   Bool_t valid;
   Int_t isRsnMini = AliAnalysisManager::GetGlobalInt("rsnUseMiniPackage",valid);
   Int_t useMC = AliAnalysisManager::GetGlobalInt("rsnUseMC",valid);
      
   if (isRsnMini) {
      postfix.Prepend("Mini");
      AliRsnMiniAnalysisTask *taskRsnMini = new AliRsnMiniAnalysisTask(Form("Rsn%s",postfix.Data()),useMC);
      Int_t refreshPrint = AliAnalysisManager::GetGlobalInt("rsnMixPrintRefresh",valid);
      if (valid) taskRsnMini->SetMixPrintRefresh(refreshPrint);
      task = (AliAnalysisTaskSE *) taskRsnMini;
   }
   else {
      AliRsnAnalysisTask *taskRsn = new AliRsnAnalysisTask(Form("Rsn%s",postfix.Data()));
      task = (AliAnalysisTaskSE *) taskRsn;
   }

   // TODO this is tmp hack
   if (!rsnIH) rsnIH = new AliRsnInputHandler();

   //   gROOT->LoadMacro("RsnConfig.C");
   if (!RsnConfig(task,useMC,isMixing,rsnIH,listRsn)) {
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
   mgr->ConnectOutput(task, 1, output);

   return task;
}
