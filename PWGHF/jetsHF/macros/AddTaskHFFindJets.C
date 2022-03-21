AliAnalysisTaskHFFindJets *AddTaskHFFindJets(TString suffix="",
                                             TString jsonconfigfile="",
                                             Bool_t readMC=kFALSE)
{
  // Creates, configures and attaches to the train the task for tracking checks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFFindJets", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFFindJets", "This task requires an input event handler");
    return NULL;
  }
  AliInputEventHandler* h= (AliInputEventHandler*)mgr->GetInputEventHandler();
  h->SetNeedField(kTRUE);

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskHFFindJets", "This task requires to run on ESD");
    return NULL;
  }
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    if (!handler) {
      ::Error("AddTaskHFFindJets","Macro called with readMC=true but MC handler not present");
      return 0;
    }
  }
  
  // Create and configure the task
  AliAnalysisTaskHFFindJets *tasktr = new AliAnalysisTaskHFFindJets();
  if(jsonconfigfile.Contains("alien://")){
    Bool_t ok=TFile::Cp(jsonconfigfile.Data(),"local_json.txt");
    if(!ok){
      ::Error("AddTaskHFFindJets","Copy of JSON file from alien failed");
      jsonconfigfile="";
    }else jsonconfigfile="local_json.txt";
  }
  if (jsonconfigfile != "") tasktr->InitFromJson(jsonconfigfile);
  tasktr->SetReadMC(readMC);

  // Add your task to the manager
  mgr->AddTask(tasktr);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  // Resolve the name of the output file
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":HFFindJets"; // create a subfolder in the file

  TString listname="clistHFFindJets";
  listname+=suffix.Data();

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(listname,
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            outputFileName.Data() );

  
  // Connect the manager to your task
  mgr->ConnectInput(tasktr,0, mgr->GetCommonInputContainer());
  
  // Same for the output
  mgr->ConnectOutput(tasktr, 1, coutput);

  // Return a pointer to the task
  return tasktr;
}
