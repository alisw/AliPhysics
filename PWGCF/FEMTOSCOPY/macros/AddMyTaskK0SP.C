//used to instantiate an object of the task,define input and output and connect it to manager

AliAnalysisTaskK0SPFemto* AddMyTaskK0SP(TString name = "name"){

  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the fun2ction
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddMyTaskK0SP", "No analysis manager found.");
    return 0x0;
  }
  // get the input event handler, again via a static method. 
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddMyTaskK0SP", "This task requires an input event handler");
    return 0x0;
  }

  // by default, a file is open for writing. here, we get the filename
  TString fileName;  
  fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MyTask";      // create a subfolder in the file

  // now we create an instance of your task
  //  cout<<"--------------------------------Add task"<<endl;
  //  AliAnalysisTaskK0SPFemto *task = new AliAnalysisTaskK0SPFemto();
  AliAnalysisTaskK0SPFemto* task = new AliAnalysisTaskK0SPFemto(name.Data());   
  //  cout<<"-------------------------------- task "<<task<<endl;
  if(!task){
    Error("AddTaskK0SPFemto","AliAnalysisTaskK0SPFemto not created!");
    return 0x0;
  }
  cout<<"DENTRO addmytask------------------------------------- 1"<<endl;
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  cout<<"DENTRO addmytask------------------------------------- 2"<<endl;
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer1", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer("MyOutputContainer2", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,3,mgr->CreateContainer("MyOutputContainer3", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}
