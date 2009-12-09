AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(char* bRec = "jets",char* bGen = "jetsAODMC_UA104",UInt_t filterMask = 16)
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
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer("pwg4spec2", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_spec2_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen));

   mgr->ConnectInput  (pwg4spec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4spec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );
   
   return pwg4spec;
}
