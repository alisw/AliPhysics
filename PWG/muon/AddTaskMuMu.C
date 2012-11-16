///
/// Configure a task to get invariant mass spectrum of dimuons
///
/// author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMu(const char* outputname, TList* triggerClassesToConsider, const char* beamYear, TArrayF* centralities)                                    
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
  
  AliAnalysisTaskMuMu* task;
  
  if ( triggerClassesToConsider ) 
  {
    task = new AliAnalysisTaskMuMu((inputDataType=="ESD"),triggerClassesToConsider,beamYear,centralities);    
  }
  else 
  {
    task = new AliAnalysisTaskMuMu((inputDataType=="ESD"),beamYear,centralities);    
  }
  
  task->AddEventCut("ALL",AliAnalysisTaskMuMu::kEventAll);
  
//  task->AddEventCut("PSALL",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS);  
//  task->AddEventCut("PSALLZSPD",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventZSPD);

  task->AddEventCut("PSALLV0ANDZSPD",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventV0AND | AliAnalysisTaskMuMu::kEventZSPD);  

  task->AddEventCut("PSALLTVXZSPD",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventTVX | AliAnalysisTaskMuMu::kEventZSPD);

  task->AddEventCut("PSALLTVXV0ANDZSPD",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventTVX | AliAnalysisTaskMuMu::kEventV0AND | AliAnalysisTaskMuMu::kEventZSPD);

  
  task->AddSingleCut("MATCHLOWRABSDCA",
                     AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kDCA);
  task->AddSingleCut("MATCHLOWRABSDCAP10",
                     AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kDCA|AliAnalysisTaskMuMu::kP10);

  task->AddPairCut("MATCHLOWRABSDCABOTH",
                   AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kDCA,
                   AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kDCA);

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

