AliAnalysisTaskHFSimpleVertices *AddTaskHFSimpleVertices(TString suffix="",
							 TString jsonconfigfile="",
							 Bool_t readMC=kFALSE)
{

  // Creates, configures and attaches to the train the task for tracking checks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFSimpleVertices", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFSimpleVertices", "This task requires an input event handler");
    return NULL;
  }   
  AliInputEventHandler* h= (AliInputEventHandler*)mgr->GetInputEventHandler();
  h->SetNeedField(kTRUE);

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskHFSimpleVertices", "This task requires to run on ESD");
    return NULL;
  }
  
  //Bool_t isMC=kFALSE;
  //if (mgr->GetMCtruthEventHandler()) isMC=kTRUE;
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    if (!handler) {
      ::Error("AddTaskHFSimpleVertices","Macro called with readMC=true but MC handler not present");
      return 0;
    }
  }
  // Create and configure the task
  AliAnalysisTaskHFSimpleVertices *tasktr = new AliAnalysisTaskHFSimpleVertices();
  if(jsonconfigfile.Contains("alien://")){
    Bool_t ok=TFile::Cp(jsonconfigfile.Data(),"local_json.txt");
    if(!ok){
      ::Error("AddTaskHFSimpleVertices","Copy of JSON file from alien failed");
      jsonconfigfile="";
    }else jsonconfigfile="local_json.txt";
  }
  if (jsonconfigfile != "") tasktr->InitFromJson(jsonconfigfile);

  
  mgr->AddTask(tasktr);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":HFVertices";

  TString listname="clistHFVertices";
  listname+=suffix.Data();

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName.Data() );
  
  mgr->ConnectInput(tasktr, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(tasktr, 1, coutput1);
  
  return tasktr;
}   

