void AddTaskCPQA(Bool_t useMC = kFALSE, Bool_t useTender = kFALSE) {

   Printf("adding task QA for central production\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskCPQA", "No analysis manager to connect to.");
      return NULL;
   }

   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskCPQA", "This task requires an input event handler.");
      return NULL;
   }

   
  
  AliAnalysisTaskCPQA *task = new AliAnalysisTaskCPQA("QADiffractive");
  task->UseMC(useMC);
   

if(useTender) {
  gROOT->LoadMacro("$ALICE_ROOT/TENDER/TenderSupplies/AddTaskTender.C");
  gROOT->LoadMacro("AddTaskTender.C");
  AliAnalysisTask* tender=0x0;
if(!useMC)
    {
        tender = AddTaskTender(kTRUE);
      // tender->SetDebugLevel(10);
    }
  else
    {
      tender = AddTaskTender(kFALSE);
      // tender->SetDebugLevel(10);
    }
 }


 if(!useMC && 0)
   {
    
     gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");

    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(useMC);
    physSelTask->GetPhysicsSelection()->SetUseBXNumbers(kFALSE);

    if(!useMC){
      
    gROOT->LoadMacro("ConfigurePhysSelection.C");
    AliOADBPhysicsSelection * oadb =  ConfigurePhysSelection(runNb);

    AliOADBFillingScheme * fsDefault = 0;

    physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadb,fsDefault);
  
      // task->SelectCollisionCandidates(AliVEvent::kUserDefined);
    }
       task->SelectCollisionCandidates(AliVEvent::kMB);
   }

 TString basefilename = AliAnalysisManager::GetCommonFileName();

   AliAnalysisDataContainer *output = mgr->CreateContainer("UDQAClist", TList::Class(), AliAnalysisManager::kOutputContainer, basefilename.Data());

  // add our task to the manager
  mgr->AddTask(task);

  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);


   return;
 }
