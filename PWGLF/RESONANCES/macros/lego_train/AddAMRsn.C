#ifndef __CINT__
#include "AliRsnTrainManager.h"
#endif
Bool_t AddAMRsnTrain(TString analysisSource = "proof", TString analysisMode = "test",TString input="aod",TString inputMC="", TString postfix = "",TString idStr="0")
{

   Bool_t usePrivateTrain = kFALSE;
   usePrivateTrain = kTRUE;

   TString legoTrainPath = "$ALICE_PHYSICS/PWGLF/RESONANCES/macros/lego_train";
//   legoTrainPath = "/home/mvala/git/AliRsn/PWGLF/RESONANCES/macros/lego_train";
   AliAnalysisManager::SetGlobalStr("RsnLegoTrainPath",legoTrainPath.Data());

   AliAnalysisManager *mrg = AliAnalysisManager::GetAnalysisManager();

   TString rsnBaseSettings = "Rsn_pp";
//   rsnBaseSettings = "Rsn_PbPb";
//   rsnBaseSettings = "Rsn_pPb";

   Bool_t useRsnMini = kTRUE;
//   useRsnMini = kFALSE;


   Bool_t useMixing = kFALSE;
//   useMixing = kTRUE;

   // RSN Setting (same as old AddRsnToManager<Rsn>.C)
   // Rsn Particle
   TString rsnStr="Phi";
   // Rsn Cut
   TString rsnCutStr="";

   rsnCutStr="PhiNsigma:KTPCnsig30";

   if ((rsnCutStr.IsNull())&&(!postfix.IsNull())) {
      rsnCutStr = "PhiNsigma:";
      rsnCutStr.Append(postfix.Data());
   }

   // Rsn Quality Cut
   TString rsnQualityCutStr = "";
//   rsnQualityCutStr = "pp_LHC11a_p4_120";
//   rsnQualityCutStr = "pp_LHC11a_p4_70";


   TString extraMacro = "";
   TString extraMacroArgs = "";
//   extraMacro = "RsnTrainSettingsExtra.C";
//   extraMacroArgs = "10.0,10,1,1,1,1,1,1,1,0";
//   extraMacroArgs = "10, 5, 5, -1, 1, 0, 1, 1, 1, 0";

   input.ToLower();
   inputMC.ToLower();
   Bool_t useMC = !inputMC.CompareTo("mc");
   Bool_t valid;

   if (usePrivateTrain) {
      if (!RsnLoadMacro("RsnPrivateTrainBaseSettings.C")) return kFALSE;

      RsnPrivateTrainBaseSettings();

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
      Int_t useVZEROEPSelection = AliAnalysisManager::GetGlobalInt("rsnUseVZEROEPSelection",valid);
      Int_t usePIDqa = AliAnalysisManager::GetGlobalInt("rsnUsePIDqa",valid);

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
            Printf("Adding PIDResponse task ...");
            gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
            AddTaskPIDResponse(useMC);
         }
      }

      if (multiInputHandler && useRsnIH) {
         // add Rsn input handler (it has to be after ESD,MC,Tender input handler, but before Mixing)
         rsnIH = new AliRsnInputHandler();
         multiInputHandler->AddInputEventHandler(rsnIH);
      }

      if (physSel>0) {
         if (!input.CompareTo("esd")) {
            gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
            Bool_t physSelBigOut = kTRUE;
//            physSelBigOut = kFALSE;

            AddTaskPhysicsSelection(useMC,kTRUE,0,physSelBigOut);
            if (physSelBigOut) mrg->SetSpecialOutputLocation("root://aaa//aaa/");
         }

         // maybe we can put it in $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C
         if (multiInputHandler) {
            AliInputEventHandler *ih = multiInputHandler->GetFirstInputEventHandler();
            ih->SetEventSelection(multiInputHandler->GetEventSelection());
         }
      }

      if (useCentralityTask) {
         gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
         AliCentralitySelectionTask *centralityTask = AddTaskCentrality(kFALSE);
      }

      if (useEventPlaneTask) {
         gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
         AliEPSelectionTask *eventPlaneTask = AddTaskEventplane();
      }

      if (useVZEROEPSelection) {
         gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
         AddTaskVZEROEPSelection();
      }

      if (usePIDqa) {
         gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
         AddTaskPIDqa();
      }

   } else {

      gSystem->Load("libEventMixing");
      gSystem->Load("libCORRFW");
      gSystem->Load("libPWGLFresonances");
   }

   if (!input.CompareTo("esd")) rsnBaseSettings.Append("_ESD");
   else rsnBaseSettings.Append("_AOD");

   // use mc
   if (useMC) rsnBaseSettings.Append("_MC");

   // use mini
   if (useRsnMini) rsnBaseSettings.Append("_MINI");

   // use mixing
   if (useMixing) rsnBaseSettings.Append("_MIX");

   if (!RsnLoadMacro("AddRsnTaskTrain.C")) return kFALSE;
   AddRsnTaskTrain(rsnBaseSettings.Data(),rsnStr.Data(),rsnCutStr.Data(),rsnQualityCutStr.Data(),extraMacro,extraMacroArgs);

   Printf("%s_%s_%s %s",rsnBaseSettings.Data(),rsnStr.Data(),rsnCutStr.Data(),rsnQualityCutStr.Data());

   return kTRUE;
}

Bool_t RsnLoadMacro(TString macro,TString path="") {

   Bool_t valid;
   TString lego_path = AliAnalysisManager::GetGlobalStr("RsnLegoTrainPath",valid);
   if (!valid) lego_path = "$ALICE_PHYSICS/PWGLF/RESONANCES/macros/lego_train";

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
