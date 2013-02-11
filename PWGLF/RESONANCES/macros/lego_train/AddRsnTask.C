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
   Int_t isMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);
   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t isAOD049 = AliRsnTrainManager::GetGlobalInt("RsnUseAOD049Patch",valid);

   if (isRsnMini) {
      postfix.Prepend("Mini");
      AliRsnMiniAnalysisTask *taskRsnMini = new AliRsnMiniAnalysisTask(TString::Format("Rsn%s",postfix.Data()).Data(),useMC);
      Int_t refreshPrint = AliRsnTrainManager::GetGlobalInt("RsnMixPrintRefresh",valid);
      if (valid) taskRsnMini->SetMixPrintRefresh(refreshPrint);
      task = (AliAnalysisTaskSE *) taskRsnMini;
   }
   else {
      AliRsnAnalysisTask *taskRsn = new AliRsnAnalysisTask(TString::Format("Rsn%s",postfix.Data()).Data());
      task = (AliAnalysisTaskSE *) taskRsn;
   }

   postfix.Append(TString::Format("_%s_%s",rsnPart.Data(),rsnCut.Data()).Data());

   if (physSelBit>=0) task->SelectCollisionCandidates((AliVEvent::EOfflineTriggerTypes)physSelBit);

   AliRsnInputHandler *rsnIH=0;
   // TODO this is tmp hack
   if (!rsnIH) rsnIH = new AliRsnInputHandler();


   TList *listRsn = new TList();
   listRsn->Add(new TNamed(rsnPart.Data(),rsnCut.Data()));

   if (!RsnLoadMacroTask("RsnConfig.C")) return 0;
   if (!RsnConfig(task,rsnIH,listRsn)) {
      Printf("Error in RsnConfig.C");
      return 0;
   }

   // setup Event Mixing
   if (isMixing) AddEventMixingSettings(task);

   if (isAOD049 && (!useMC) && (collisionType==1)) {
      if (isRsnMini) taskRsnMini->SetUseCentralityPatch(kTRUE);
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

   //just small fix to replace ':' by '_' (easier to read)
   postfix.ReplaceAll(":","_");
   // create containers for output
   AliAnalysisDataContainer *output = mgr->CreateContainer(TString::Format("RsnHist%s", postfix.Data()).Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath.Data());

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

void AddEventMixingSettings(AliAnalysisTaskSE *task) {

   Bool_t valid = kTRUE;
   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t mixNum = AliRsnTrainManager::GetGlobalInt("RsnNumMix",valid);

   Double_t mixDiffMult = AliRsnTrainManager::GetGlobalInt("RsnMixDiffMult",valid);
   Double_t mixDiffVz = AliRsnTrainManager::GetGlobalInt("RsnMixDiffVz",valid);
   Double_t mixDiffAngle = AliRsnTrainManager::GetGlobalInt("RsnMixDiffAngle",valid);

   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsn = (AliRsnMiniAnalysisTask *) task;
      if (collisionType == 0) {
         //         taskRsn->UseMultiplicity("TRACKS");
         taskRsn->UseMultiplicity("QUALITY");
      } else {
         taskRsn->UseCentrality("V0M");
      }
      if (mixDiffMult>0.0) taskRsn->SetMaxDiffMult(mixDiffMult);

      // set mixing
      taskRsn->UseContinuousMix();
      //task->UseBinnedMix();
      taskRsn->SetNMix(mixNum);

      if (mixDiffVz>0.0) taskRsn->SetMaxDiffVz(mixDiffVz);
      if (mixDiffAngle>0.0) taskRsn->SetMaxDiffAngle(mixDiffAngle);
      // 30.0 * TMath::DegToRad()
   }
   // TODO RSN non Mini

}

Bool_t RsnLoadMacroTask(TString macro,TString path="") {

   Bool_t valid;
   TString lego_path = AliAnalysisManager::GetGlobalStr("RsnLegoTrainPath",valid);
   if (!valid) lego_path = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";

   if (!gSystem->AccessPathName(macro.Data())) {
      gROOT->LoadMacro(macro.Data());
      Printf("Macro loaded from %s/%s ...",gSystem->pwd(),macro.Data());
      return kTRUE;
   }

   if (!gSystem->AccessPathName(gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())))) {
      gROOT->LoadMacro(gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())));
      Printf("Macro loaded from %s ...",gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())));
      return kTRUE;
   }

   Printf("Error loading %s",macro.Data());

   return kFALSE;
}