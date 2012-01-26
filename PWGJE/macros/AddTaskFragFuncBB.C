
AliAnalysisTaskFragFuncBB *AddTaskFragFuncBB(UInt_t filterMask=16, Int_t iPhysicsSelection=1){
	AddTaskFragFuncBB("jets", "", filterMask, iPhysicsSelection);
}

AliAnalysisTaskFragFuncBB *AddTaskFragFuncBB(UInt_t filterMask=16, Bool_t kUseAODMC=kFALSE, Int_t iPhysicsSelection=1, UInt_t iFlag){
    AliAnalysisTaskFragFuncBB *ff=0;
	if(kUseAODMC){
	    if(iFlag&(1<<0)) ff = AddTaskFragFuncBB("jets", "jetsAODMC2b_UA104", filterMask, iPhysicsSelection);
		if(iFlag&(1<<1)) ff = AddTaskFragFuncBB("jetsAOD_UA107", "jetsAODMC2b_UA107", filterMask, iPhysicsSelection);
	}
	if(iFlag&(1<<2)) ff = AddTaskFragFuncBB("jets", "", filterMask, iPhysicsSelection);
	
	return ff;
}

// _______________________________________________________________________________________

AliAnalysisTaskFragFuncBB *AddTaskFragFuncBB(
	const char* bRecJets,
	const char* bGenJets,
	UInt_t filterMask,
	Int_t iPhysicsSelection)
{
   // Creates a fragmentation function task,
   // configures it and adds it to the analysis manager.
   
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
	  ::Error("AddTaskFragFuncBB", "No analysis manager to connect to.");
	  return NULL;
   }
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
	 ::Error("AddTaskFragFunc", "This task requires an input event handler");
	  return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   TString typeRec(bRecJets);
   TString typeGen(bGenJets);
   typeRec.ToUpper();
   typeGen.ToUpper();
   
   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskFragFuncBB *fragfunc = new AliAnalysisTaskFragFuncBB(Form("Fragmenation Function %s %s", bRecJets, bGenJets));
   
   // pwg4fragfunc->SetAnalysisType(AliAnalysisTaskFragFunc::kAnaMC);
   // if(iAODanalysis) pwg4fragfunc->SetAODInput(kTRUE);
   // pwg4fragfunc->SetDebugLevel(11);
   
   Printf("Rec Jets %s", bRecJets);
   Printf("Gen Jets %s", bGenJets);
   
   fragfunc->SetBranchGenJets(bGenJets);
   fragfunc->SetBranchRecJets(bRecJets);
   
   fragfunc->SetFilterMask(filterMask);
   //fragfunc->SetUseGlobalSelection(kTRUE);
   
   if(type == "AOD"){
	  // Assume all jets are produced already
	  fragfunc->SetAODJetInput(kTRUE);
	  fragfunc->SetAODTrackInput(kTRUE);
	  fragfunc->SetAODMCInput(kTRUE);
   }
	
   if(typeRec.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
	 fragfunc->SetTrackTypeRec(AliAnalysisTaskFragFunc::kTrackAODMCChargedAcceptance);
   }
   else if (typeRec.Contains("AODMC2")){
	 fragfunc->SetTrackTypeRec(AliAnalysisTaskFragFunc::kTrackAODMCCharged);
   }
   else if (typeRec.Contains("AODMC")){
	 fragfunc->SetTrackTypeRec(AliAnalysisTaskFragFunc::kTrackAODMCAll);
   }
   else { // catch akk use AOD
	 fragfunc->SetTrackTypeRec(AliAnalysisTaskFragFunc::kTrackAOD);
   }

   if(typeGen.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
	 fragfunc->SetTrackTypeGen(AliAnalysisTaskFragFunc::kTrackAODMCChargedAcceptance);
   }
   else if (typeGen.Contains("AODMC2")){
	 fragfunc->SetTrackTypeGen(AliAnalysisTaskFragFunc::kTrackAODMCCharged);
   }
   else if (typeGen.Contains("AODMC")){
	 fragfunc->SetTrackTypeGen(AliAnalysisTaskFragFunc::kTrackAODMCAll);
   }
   else if (typeGen.Length()>0){ // catch all use AOD
	 fragfunc->SetTrackTypeGen(AliAnalysisTaskFragFunc::kTrackAOD);
   }
   else {  //
	  fragfunc->SetTrackTypeGen(AliAnalysisTaskFragFunc::kTrackKineCharged);
   }

   //if(iPhysicsSelection)fragfunc->SelectCollisionCandidates();

   mgr->AddTask(fragfunc);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_FragFunc = mgr->CreateContainer(Form("fracfuncBB_%s_%s",bRecJets,bGenJets), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_fragfuncBB_%s_%s",AliAnalysisManager::GetCommonFileName(),bRecJets,bGenJets));

   mgr->ConnectInput  (fragfunc, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (fragfunc, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (fragfunc, 1, coutput1_FragFunc);
   
   return fragfunc;
}