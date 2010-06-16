AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(const char* bRec = "jets",const char* bGen = "jetsAODMC_UA104",UInt_t filterMask = 16, Int_t iPhysicsSelection = 1);


AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2Delta(UInt_t filterMask = 16,Bool_t kUseAODMC = kFALSE,Int_t iPhysicsSelection = 1,UInt_t iFlag){
  AliAnalysisTaskJetSpectrum2 *js = 0;
  if(kUseAODMC){
    if(iFlag&(1<<0))js = AddTaskJetSpectrum2("jets","jetsAODMC_UA104",filterMask,iPhysicsSelection);
    if(iFlag&(1<<1))js = AddTaskJetSpectrum2("jets","jetsAODMC2_UA104",filterMask,iPhysicsSelection);

    if(iFlag&(1<<2))js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","jetsAODMC_FASTJET04",filterMask,iPhysicsSelection);
    if(iFlag&(1<<3))js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","jetsAODMC2_FASTJET04",filterMask,iPhysicsSelection);

    if(iFlag&(1<<4)){
      js = AddTaskJetSpectrum2("jetsAOD_FASTKT04","jetsAODMC_FASTKT04",filterMask,iPhysicsSelection);
    }
    if(iFlag&(1<<5))js = AddTaskJetSpectrum2("jetsAOD_FASTKT04","jetsAODMC2_FASTKT04",filterMask,iPhysicsSelection);
    if(iFlag&(1<<6))js = AddTaskJetSpectrum2("jetsAOD_UA107","jetsAODMC_UA107",filterMask,iPhysicsSelection);
  }

  if(iFlag&(1<<7))js = AddTaskJetSpectrum2("jets","jetsAOD_FASTJET04",filterMask,iPhysicsSelection);

  if(iFlag&(1<<8))js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","",filterMask,iPhysicsSelection);
  if(iFlag&(1<<9))js = AddTaskJetSpectrum2("jetsAOD_FASTKT04","",filterMask,iPhysicsSelection);
  if(iFlag&(1<<10))js = AddTaskJetSpectrum2("jetsAOD_SISCONE04","",filterMask,iPhysicsSelection);
  
  if(iFlag&(1<<11)){
    js = AddTaskJetSpectrum2("jetsAOD_UA107","",filterMask,iPhysicsSelection);
    js->SetRecEtaWindow(0.2);
  }
  if(iFlag&(1<<12)){
    js = AddTaskJetSpectrum2("jetsAOD_FASTJET07","",filterMask,iPhysicsSelection);
    js->SetRecEtaWindow(0.2);
  }
  if(iFlag&(1<<13)){
    js = AddTaskJetSpectrum2("jetsAOD_FASTKT07","",filterMask,iPhysicsSelection);
    js->SetRecEtaWindow(0.2);
  }
  if(iFlag&(1<<14)){
    js = AddTaskJetSpectrum2("jetsAOD_SISCONE07","",filterMask,iPhysicsSelection);
    js->SetRecEtaWindow(0.2);
  }
  return js;
}


AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(const char* bRec,const char* bGen ,UInt_t filterMask,Int_t iPhysicsSelection)
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
   TString typeRec(bRec);
   TString typeGen(bGen);
   typeGen.ToUpper();
   typeRec.ToUpper();
   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskJetSpectrum2* pwg4spec = new  AliAnalysisTaskJetSpectrum2(Form("JetSpectrum2%s%s",bRec,bGen));
      
   // or a config file
   // pwg4spec->SetAnalysisType(AliAnalysisTaskJetSpectrum2::kAnaMC);
   // if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   // pwg4spec->SetDebugLevel(11); 

   pwg4spec->SetBranchGen(bGen); 
   pwg4spec->SetBranchRec(bRec); 


   pwg4spec->SetFilterMask(filterMask); 
   pwg4spec->SetUseGlobalSelection(kTRUE); 
   pwg4spec->SetMinJetPt(5.);


   if(type == "AOD"){
     // Assume all jet are produced already
     pwg4spec->SetAODJetInput(kTRUE);
     pwg4spec->SetAODTrackInput(kTRUE);
     pwg4spec->SetAODMCInput(kTRUE);
   }

   if(typeRec.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODMCChargedAcceptance);
   }
   else if (typeRec.Contains("AODMC2")){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODMCCharged);
   }
   else if (typeRec.Contains("AODMC")){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODMCAll);
   }
   else { // catch akk use AOD
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAOD);
   }

   if(typeGen.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCChargedAcceptance);
   }
   else if (typeGen.Contains("AODMC2")){
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCCharged);
   }
   else if (typeGen.Contains("AODMC")){
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCAll);
   }
   else if (typeGen.Length()>0){ // catch all use AOD
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAOD);
   }

   if(iPhysicsSelection)pwg4spec->SelectCollisionCandidates();

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
