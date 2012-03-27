Bool_t AddAMEventMixingTest(TString analysisSource = "proof", TString analysisMode = "test",TString input="aod",TString inputMC="", TString postfix = "",TString idStr="0")
{
  
   Bool_t useEventMixingPar      = 0;

   Int_t usePhysSel              = 0;
   
   Bool_t useMC = !inputMC.CompareTo("mc");

   // ALICE stuff
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error[AddAMEventMixingTest] mgr is null !!!"); return kFALSE; }
   
   
   AliAnalysisGrid *analysisPlugin = mgr->GetGridHandler();
   if (!analysisPlugin) { Printf("Error[AddAMEventMixingTest] : analysisPlugin is null !!!"); return kFALSE; }

   TString myAdditionalLibs;
   if (useEventMixingPar) { AliAnalysisAlien::SetupPar("EventMixing"); myAdditionalLibs += " EventMixing.par"; }
   else { gSystem->Load("libEventMixing.so"); myAdditionalLibs += " libEventMixing.so"; }
   
   gROOT->LoadMacro("AliAnalysisTaskEx02.cxx++g");
   analysisPlugin->SetAnalysisSource("AliAnalysisTaskEx02.cxx+");
   myAdditionalLibs+=" AliAnalysisTaskEx02.h AliAnalysisTaskEx02.cxx";
   analysisPlugin->SetAdditionalLibs(myAdditionalLibs.Data());
   

  AliMultiInputEventHandler *multiInputHandler = mgr->GetInputEventHandler();

   if (usePhysSel) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AddTaskPhysicsSelection(useMC);

      // maybe we can put it in $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C
      AliMultiInputEventHandler *multiIH = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
      if (multiIH){
         AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler *>(multiIH->GetFirstInputEventHandler());
         if (esdIH) esdIH->SetEventSelection(multiIH->GetEventSelection());
         AliAODInputHandler *aodIH = dynamic_cast<AliAODInputHandler *>(multiIH->GetFirstInputEventHandler());
         if (aodIH) aodIH->SetEventSelection(multiIH->GetEventSelection());
      }
   }


  // add mixing handler (uncomment to turn on Mixnig)
  gROOT->LoadMacro("AddMixingHandler.C");
  AddMixingHandler(multiInputHandler, input, useMC,postfix);
  
   

//    // load and run AddTask macro
   gROOT->LoadMacro("AddEventMixingTestTask.C");
   AddEventMixingTestTask(input, useMC, postfix);

   return kTRUE;
}
