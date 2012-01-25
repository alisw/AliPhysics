AliAnalysisTaskJetBackgroundSubtract *AddTaskJetBackgroundSubtract(TString sJetBranches,Int_t iSubtraction = 1,char* cToReplace = "B0",char* cReplaceMask = "B%d", char* cSuffix = ""){




  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_obusch_jets", "No analysis manager found.");
    return 0;
  }




  //========= Add task to the ANALYSIS manager =====
  
  // TOkenize...
  TObjArray *objArr = sJetBranches.Tokenize(" ");


  TString cAdd;
  cAdd += Form(cReplaceMask,iSubtraction);
  if(cSuffix!="") cAdd += Form("_%s",cSuffix);

  AliAnalysisTaskJetBackgroundSubtract *task = new AliAnalysisTaskJetBackgroundSubtract(Form("JetSubtract_%s",cAdd.Data()));
  for(int iJB = 0;iJB<objArr->GetEntries();iJB++){
    TObjString *ostr = (TObjString*)objArr->At(iJB);
    task->AddJetBranch(ostr->GetString().Data());
  }
  task->SetToReplace(cToReplace); 
  task->SetReplacementMask(cReplaceMask); 
  task->SetSubtractionMethod(iSubtraction);
  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();	
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(
							   Form("pwgJetSubtract_%s",cAdd.Data()),
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer,
							   Form("%s:PWG4_JetSubtract_%s",AliAnalysisManager::GetCommonFileName(),cAdd.Data()));

  mgr->ConnectInput(task,0,cinput );
  mgr->ConnectOutput(task,1,coutput);
 
  return task;
}
