#ifndef __CINT__
#include <AliAnalysisManager.h>
#include <AliLog.h>
#endif

Bool_t AddAMRsn(TString analysisSource = "proof", TString analysisMode = "test",TString input="aod",TString inputMC="", TString postfix = "",TString idStr="0")
{

   analysisSource.ToLower(); analysisMode.ToLower();

   if (!RsnLoadMacro("RsnManager.C")) return kFALSE;
   TList *listRsn = RsnManager();

   Bool_t useMC = !inputMC.CompareTo("mc");
   input.ToLower();
   Bool_t valid;

   Int_t eventMixinPar = AliAnalysisManager::GetGlobalInt("rsnUseEventMixingPar",valid);
   Int_t rsnPar = AliAnalysisManager::GetGlobalInt("rsnUseRSNPar",valid);
   Int_t rsnParDev = AliAnalysisManager::GetGlobalInt("rsnUseRSNParDev",valid);
   if (eventMixinPar) rsnPar = 1;
   if (rsnPar&&rsnParDev>=0) rsnParDev=1;

   Int_t pidResponse = AliAnalysisManager::GetGlobalInt("rsnUsePIDResponse",valid);
   Int_t useRsnIH = AliAnalysisManager::GetGlobalInt("rsnUseRsnInputHandler",valid);
   Int_t physSel = AliAnalysisManager::GetGlobalInt("rsnUsePhysSel",valid);
   Int_t useCentralityTask = AliAnalysisManager::GetGlobalInt("rsnUseCentralityTask",valid);
   Int_t useEventPlaneTask = AliAnalysisManager::GetGlobalInt("rsnUseEventPlaneTask",valid);
   Int_t usePIDqa = AliAnalysisManager::GetGlobalInt("rsnUsePIDqa",valid);
   Int_t splitMgrByTask = AliAnalysisManager::GetGlobalInt("rsnSplitMgrByTasks",valid);

   Int_t useMixing = AliAnalysisManager::GetGlobalInt("rsnUseMixing",valid);

   Int_t isRsnMini = AliAnalysisManager::GetGlobalInt("rsnUseMiniPackage",valid);
   Int_t mixNum = AliAnalysisManager::GetGlobalInt("rsnNumMix",valid);

   // ALICE stuff
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error[AddAMRsn] mgr is null !!!"); return kFALSE; }

   AliAnalysisGrid *analysisPlugin = mgr->GetGridHandler();
   if (!analysisPlugin) { Printf("Error[AddAMRsn] : analysisPlugin is null !!!"); return kFALSE; }

   TString myAdditionalLibs;
   if (eventMixinPar) { AliAnalysisAlien::SetupPar("EventMixing"); myAdditionalLibs += " EventMixing.par"; }
   else { gSystem->Load("libEventMixing"); myAdditionalLibs += " libEventMixing.so"; }

   TString rsnLibName = "PWGLFresonances";
   if (gSystem->Getenv("ALICE_ROOT")) {
      TString alirootVersion = gSystem->GetFromPipe("aliroot --version | awk '{print $3}'");
      if (alirootVersion<"v5-02-19-AN" && alirootVersion.CompareTo("trunk")) rsnLibName = "PWG2resonances";
      if (rsnPar) { AliAnalysisAlien::SetupPar(rsnLibName.Data()); myAdditionalLibs += Form(" %s.par",rsnLibName.Data()); }
      else { gSystem->Load(Form("lib%s",rsnLibName.Data())); myAdditionalLibs += Form(" lib%s.so",rsnLibName.Data()); }
   }

   if (rsnParDev>=0) {
      if (rsnParDev) { AliAnalysisAlien::SetupPar("PWGLFresonancesdev"); myAdditionalLibs += " PWGLFresonancesdev.par"; }
      else { gSystem->Load("libPWGLFresonancesdev"); myAdditionalLibs += " libPWGLFresonancesdev.so"; }
   }
   analysisPlugin->SetAdditionalLibs(myAdditionalLibs.Data());

   AliMultiInputEventHandler *multiInputHandler = 0;
   AliInputEventHandler *inputHandler = mgr->GetInputEventHandler();

   TString className = inputHandler->ClassName();
   if (!className.CompareTo("AliMultiInputEventHandler")) {
      multiInputHandler = (AliMultiInputEventHandler *)inputHandler;
   }

   AliRsnInputHandler *rsnIH=0;

   if (pidResponse) {
      if (multiInputHandler) {
         // add PID Response Handler
         if (!RsnLoadMacro("AddPIDResponseInputHandler.C")) return kFALSE;
         AddPIDResponseInputHandler(multiInputHandler,useMC);
      } else {
         gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
         AddTaskPIDResponse(useMC);
      }
   }

   if (multiInputHandler && useRsnIH) {
      // add Rsn input handler (it has to be after ESD,MC,Tender input handler, but before Mixing)
      rsnIH = new AliRsnInputHandler();
      multiInputHandler->AddInputEventHandler(rsnIH);
   }

   if (physSel) {
      if (!input.CompareTo("esd")) {
         gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
         AddTaskPhysicsSelection(useMC);
      }

      // maybe we can put it in $ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C
      if (multiInputHandler) {
         AliInputEventHandler *ih = multiInputHandler->GetFirstInputEventHandler();
         ih->SetEventSelection(multiInputHandler->GetEventSelection());
      }
   }

   if (useCentralityTask) {
      gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *centralityTask = AddTaskCentrality(kFALSE);
   }

   if (useEventPlaneTask) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
      AliEPSelectionTask *eventPlaneTask = AddTaskEventplane();
   }

   if (usePIDqa) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
      AddTaskPIDqa();
   }

   // load and run AddTask macro
   if (!RsnLoadMacro("AddRsnAnalysisTask.C")) return kFALSE;
   if (!RsnLoadMacro("RsnConfig.C")) return kFALSE;
   if (!RsnLoadMacro("AddMixingHandler.C")) return kFALSE;
   if (!analysisSource.CompareTo("grid")) {
      if (!RsnLoadMacro("RsnGridPlugin.C")) return kFALSE;
      RsnGridPlugin(analysisMode);
   }

   if (splitMgrByTask) {
      Int_t iTask=0;
      TList *l=0;
      TNamed *rsnObj = 0;
      AliAnalysisTaskSE *task=0;
      TString rsnName,rsnCutName;
      TIter next(listRsn);
      while ((rsnObj = (TNamed *)next())) {
         l = new TList();
         Printf("Adding task for RSN:%s CUT:%s ",rsnObj->GetName(),rsnObj->GetTitle());
         l->Add(new TNamed(*rsnObj));
         task = AddRsnAnalysisTask(input, useMC, useMixing,rsnIH,l,Form("%s_%d",postfix.Data(),iTask++));
         if (useMixing) {
            // add mixing handler (uncomment to turn on Mixnig)
            AddMixingHandler(multiInputHandler,task, input, useMC,isRsnMini, mixNum,postfix);
         }
      }
   } else {
      task = AddRsnAnalysisTask(input, useMC, useMixing,rsnIH,listRsn,postfix);
      if (useMixing) {
         // add mixing handler (uncomment to turn on Mixnig)
         AddMixingHandler(multiInputHandler,task, input, useMC,isRsnMini, mixNum,postfix);
      }
   }

   //    mgr->AddClassDebug("AliRsnCutTrackQuality",AliLog::kDebug+3);

   return kTRUE;
}

Bool_t RsnLoadMacro(TString macro,TString path="") {

   Bool_t valid;
   TString lego_path = AliAnalysisManager::GetGlobalStr("rsnLegoTrainPath",valid);
   if (!valid) lego_path = "$ALICE_ROOT/PWG2/RESONANCES/macros/lego_train";

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
