AliAnalysisMultPt* AddMultPt(const char *suffix = "")
{
  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

  // resolve the name of the output file

  TString name;
  name.Form("MulttpT%s", suffix);
    // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MultPtt";      // create a subfolder in the file

  // now we create an instance of your task
  AliAnalysisMultPt* task = new AliAnalysisMultPt(name.Data());
  if(!task) return 0x0;
  task->SelectCollisionCandidates(AliVEvent::kINT7); //to check which triggers are selected
  task->SetMCRead(kTRUE);
  task->SetPileUpRead(kFALSE);
  task->SetChi2DoF(4);
  task->SetPtLimits(0.2, 5.0);
  task->SetEtaMinLimit(-0.8);
  task->SetEtaMaxLimit(0.8);
  
  printf("Container name is %s\n",name.Data());
    
  // add your task to the manager
  mgr->AddTask(task);

  // connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer(name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

  // important: return a pointer to your task
  return task;
}
