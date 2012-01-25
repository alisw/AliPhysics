AliAnalysisTaskJetChem *AddTaskJetChem(){
        
  AliAnalysisTaskJetChem *ff=0;
   
  ff = AddTask("clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00",0,AliAnalysisTaskJetChem::kOnFly);
    
  return ff;
}

// _______________________________________________________________________________________

AliAnalysisTaskJetChem *AddTask(const char* recJetsBranch,Int_t eventClass, Int_t K0type)
{
  // Creates a jet chem task,
  // configures it and adds it to the analysis manager.
  
   
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetChem", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetChem", "This task requires an input event handler");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Printf("Data Type: %s", type.Data());
  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskJetChem *task = new AliAnalysisTaskJetChem("bla");

  Int_t debug = -1; // debug level
  if(debug>=0) task->SetDebugLevel(debug);
  
  TString branchRecJets(recJetsBranch);
  if(!branchRecJets.Contains("noRecJets")) task->SetBranchRecJets(branchRecJets);
  
  task->SetFilterMask(256); // for h+/h-
  task->SetEventClass(eventClass);
  
  task->SetK0Type(K0type);

  if(K0type == AliAnalysisTaskJetChem::kOnFlyPrim || AliAnalysisTaskJetChem::kOfflPrim) task->SetFilterMaskK0(256);

  task->SetFFRadius(); 

  task->SetTrackCuts(0.15, -0.75, 0.75, 0., 2*TMath::Pi());
  task->SetJetCuts(5., -0.35, 0.35, 0., 2*TMath::Pi());

  // Define histo bins
   task->SetFFHistoBins();
   task->SetQATrackHistoBins();
   task->SetFFInvMassHistoBins();
   task->SetPhiCorrInvMassHistoBins();
   
   mgr->AddTask(task);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   TString strK0type;
   if(K0type ==  AliAnalysisTaskJetChem::kOnFly)     strK0type = "OnFly";
   if(K0type ==  AliAnalysisTaskJetChem::kOnFlyPID)  strK0type = "OnFlyPID";
   if(K0type ==  AliAnalysisTaskJetChem::kOnFlydEdx) strK0type = "OnFlydEdx";
   if(K0type ==  AliAnalysisTaskJetChem::kOnFlyPrim) strK0type = "OnFlyPrim";
   if(K0type ==  AliAnalysisTaskJetChem::kOffl)      strK0type = "Offl";
   if(K0type ==  AliAnalysisTaskJetChem::kOfflPID)   strK0type = "OfflPID";
   if(K0type ==  AliAnalysisTaskJetChem::kOffldEdx)  strK0type = "OffldEdx";
   if(K0type ==  AliAnalysisTaskJetChem::kOfflPrim)  strK0type = "OfflPrim";
   
   TString listName(Form("PWG4_JetChem_%s_%s_cl%d",branchRecJets.Data(),strK0type.Data(),eventClass));
		    
   AliAnalysisDataContainer *coutput_FragFunc = mgr->CreateContainer(listName, 
								     TList::Class(),
								     AliAnalysisManager::kOutputContainer,
								     Form("%s:PWG4_JetChem",AliAnalysisManager::GetCommonFileName()));
   
   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, coutput_FragFunc);
   
   return task;
}
