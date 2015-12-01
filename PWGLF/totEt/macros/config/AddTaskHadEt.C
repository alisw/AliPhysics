AliAnalysisTaskHadEt* AddTaskHadEt(
   const char *outfilename       = "AnalysisOutput.root",
   Bool_t isMC                   = kFALSE,
   TString recoConfigFile       = "",
   TString mcConfigFile         = "",
   Int_t nDataSet                = 2015,
   Bool_t useDefaultConstructor  = kTRUE,
   const char *tag               = ""
)
{  
  cout<<"IN THE ADDTASK"<<endl;
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHadEt", "No analysis manager to connect to.");
    return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTaskHadEt", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name(Form("TaskHadEt%s", tag));
  // Default constructor listed below
//  AliAnalysisTaskHadEt::AliAnalysisTaskHadEt(const char *name, Bool_t isMc, TString recoConfigFile, TString mcConfigFile) :
  AliAnalysisTaskHadEt *HadEttask;
  if(useDefaultConstructor) { 
    HadEttask = new AliAnalysisTaskHadEt();
  } else {
    HadEttask = new AliAnalysisTaskHadEt(name, isMC, recoConfigFile, mcConfigFile);
  }

  // check for MC 
  if(isMC) HadEttask->SetMCData();

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(HadEttask);

  // Create containers for input/output
  mgr->ConnectInput (HadEttask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *contHadEt = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);

  // connect output container
  mgr->ConnectOutput(HadEttask,1,contHadEt);

  // set debug level
  mgr->SetDebugLevel(0);

  return HadEttask;
}
