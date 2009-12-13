AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(char* bRec = "jets",char* bGen = "jetsAODMC_UA104",UInt_t filterMask = 16);


AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2Delta(UInt_t filterMask = 16,Bool_t kUseAODMC = kFALSE){
  AliAnalysisTaskJetSpectrum2 *js = 0;
  if(kUseAODMC){
  js = AddTaskJetSpectrum2("jets","jetsAODMC_UA104",filterMask);
    js = AddTaskJetSpectrum2("jets","jetsAODMC2_UA104",filterMask);

    js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","jetsAODMC_FASTJET04",filterMask);
    js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","jetsAODMC2_FASTJET04",filterMask);

  js = AddTaskJetSpectrum2("jetsAOD_UA107","jetsAODMC_UA107",filterMask);
  }
  js = AddTaskJetSpectrum2("jets","jetsAOD_FASTJET04",filterMask);
  js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","",filterMask);
  js = AddTaskJetSpectrum2("jetsAOD_UA107","",filterMask);
  js->SetRecEtaWindow(0.2);

  return js;
}


AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(char* bRec,char* bGen ,UInt_t filterMask)
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetSpectrum2", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
     ::Error("AddTaskJetSpectrum2", "This task requires an input event handler");
      return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   TString typeMC(bGen);
   typeMC.ToUpper();
   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskJetSpectrum2* pwg4spec = new  AliAnalysisTaskJetSpectrum2(Form("Jet Spectrum %s %s",bRec,bGen));
      
   // or a config file
   // pwg4spec->SetAnalysisType(AliAnalysisTaskJetSpectrum2::kAnaMC);
   // if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   // pwg4spec->SetDebugLevel(11); 
   pwg4spec->SetBranchRec(bRec); 
   pwg4spec->SetBranchGen(bGen); 
   pwg4spec->SetFilterMask(filterMask); 
   pwg4spec->SetUseGlobalSelection(kTRUE); 

   if(type == "AOD"){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODIn);
     pwg4spec->SetAODInput(kTRUE);
   }
   else pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODOut);

   if(typeMC.Contains("AODMC2")){
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCCharged);
   }
   else if (typeMC.Contains("AODMC2")){
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCAll);
   }
   else if (typeMC.Contains("AOD")) {
     if(type == "AOD"){
       pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODIn);
     }
     else{
       pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODOut);
     }
   }
   mgr->AddTask(pwg4spec);
     
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer(Form("pwg4spec2_%s_%s",bRec,bGen), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_spec2_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen));

   mgr->ConnectInput  (pwg4spec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4spec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );
   
   return pwg4spec;
}
