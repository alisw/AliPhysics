#ifndef __CINT__
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <ANALYSIS/AliAnalysisManager.h>
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <ANALYSIS/AliAnalysisAlien.h>
#include <ANALYSIS/AliAnalysisTaskSE.h>
#endif

Bool_t AddAMRsn(TString analysisSource = "proof", TString analysisMode = "test",TString input="aod",TString inputMC="", TString postfix = "",TString idStr="0")
{

   gROOT->LoadMacro("RsnManager.C");
   TList *listRsn = RsnManager();

   Bool_t useMC = !inputMC.CompareTo("mc");

   // ALICE stuff
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error[AddAMRsn] mgr is null !!!"); return kFALSE; }

   AliAnalysisGrid *analysisPlugin = mgr->GetGridHandler();
   if (!analysisPlugin) { Printf("Error[AddAMRsn] : analysisPlugin is null !!!"); return kFALSE; }

   if (gRsnUseEventMixingPar) gRsnUseRSNPar = 1;

   TString myAdditionalLibs;
   if (gRsnUseEventMixingPar) { AliAnalysisAlien::SetupPar("EventMixing"); myAdditionalLibs += " EventMixing.par"; }
   else { gSystem->Load("libEventMixing.so"); myAdditionalLibs += " libEventMixing.so"; }

   if (gRsnUseRSNPar) { AliAnalysisAlien::SetupPar("PWG2resonances"); myAdditionalLibs += " PWG2resonances.par"; }
   else { gSystem->Load("libPWG2resonances.so"); myAdditionalLibs += " libPWG2resonances.so"; }

   analysisPlugin->SetAdditionalLibs(myAdditionalLibs.Data());


   AliMultiInputEventHandler *multiInputHandler = mgr->GetInputEventHandler();
   AliRsnInputHandler *rsnIH=0;

   if (gRsnUsePIDResponse) {
      // add PID Response Handler
      gROOT->LoadMacro("AddPIDResponseInputHandler.C");
      AddPIDResponseInputHandler(multiInputHandler);
   }

   if (gRsnUseRsnInputHandler) {
      // add Rsn input handler (it has to be after ESD,MC,Tender input handler, but before Mixing)
      AliRsnInputHandler *rsnIH = new AliRsnInputHandler();
      multiInputHandler->AddInputEventHandler(rsnIH);
   }

   if (gRsnUsePhysSel) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AddTaskPhysicsSelection(useMC);

      // maybe we can put it in $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C
      AliMultiInputEventHandler *multiIH = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
      if (multiIH) {
         AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler *>(multiIH->GetFirstInputEventHandler());
         if (esdIH) esdIH->SetEventSelection(multiIH->GetEventSelection());
         AliAODInputHandler *aodIH = dynamic_cast<AliAODInputHandler *>(multiIH->GetFirstInputEventHandler());
         if (aodIH) aodIH->SetEventSelection(multiIH->GetEventSelection());
      }
   }


   // load and run AddTask macro
   gROOT->LoadMacro("AddRsnAnalysisTask.C");
   gROOT->LoadMacro("RsnConfig.C");
   gROOT->LoadMacro("AddMixingHandler.C");

   if (gRsnSplitMgrByTasks) {
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
         task = AddRsnAnalysisTask(input, useMC, gRsnUseMixing,rsnIH,l,Form("%s_%d",postfix.Data(),iTask++));
         if (gRsnUseMixing) {
            // add mixing handler (uncomment to turn on Mixnig)
            AddMixingHandler(multiInputHandler,task, input, useMC,gRsnUseMiniPackage, gRsnNumMix,postfix);
         }
      }
   } else {
      task = AddRsnAnalysisTask(input, useMC, gRsnUseMixing,rsnIH,listRsn,postfix);
      if (gRsnUseMixing) {
         // add mixing handler (uncomment to turn on Mixnig)
         AddMixingHandler(multiInputHandler,task, input, useMC,gRsnUseMiniPackage, gRsnNumMix,postfix);
      }
   }




   return kTRUE;
}
