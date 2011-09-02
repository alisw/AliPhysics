///
/// Configure a task to get invariant mass spectrum of dimuons
///
/// author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMu(const char* outputname, TList* triggerClassesToConsider, Bool_t aa)
                                    
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Configure analysis
  //===========================================================================  
  
  AliAnalysisMuMu* task;
  
  if ( triggerClassesToConsider ) 
  {
    task = AliAnalysisMuMu::Create(inputDataType,triggerClassesToConsider);    
  }
  else 
  {
    task = AliAnalysisMuMu::Create(inputDataType,aa);    
  }
  
  task->AddSingleCut("MATCHLOWRABS",AliAnalysisMuMu::kAll|AliAnalysisMuMu::kMatchedLow|AliAnalysisMuMu::kRabs);
  task->AddPairCut("MATCHLOWRABS",AliAnalysisMuMu::kAll|AliAnalysisMuMu::kMatchedLow|AliAnalysisMuMu::kRabs);

  task->AddSingleCut("MATCHHIGHRABSDCA",AliAnalysisMuMu::kAll|AliAnalysisMuMu::kMatchedHigh|AliAnalysisMuMu::kRabs|AliAnalysisMuMu::kDCA);
  task->AddPairCut("MATCHHIGHRABSDCA",AliAnalysisMuMu::kAll|AliAnalysisMuMu::kMatchedHigh|AliAnalysisMuMu::kRabs|AliAnalysisMuMu::kDCA);

  mgr->AddTask(task);  
  
  static int n(0);
  
  ++n;
  
  TString containerName("chist");
  
  if ( n > 1 ) containerName += Form("%d",n);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();  
  AliAnalysisDataContainer *coutput = 
  mgr->CreateContainer(containerName.Data(), TList::Class(),
                       AliAnalysisManager::kOutputContainer,outputname);
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}

