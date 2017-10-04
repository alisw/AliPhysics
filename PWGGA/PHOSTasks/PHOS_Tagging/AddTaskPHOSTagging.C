AliAnalysisTaskTaggedPhotons* AddTaskPHOSTagging (const char* name = "PHOSTagging",
					    const char* options = "",
					    UInt_t offlineTriggerMask = AliVEvent::kCentral,
					    Float_t timeCut = 100.e-9, //accept clusters with |t|<timeCut
                                            Float_t distCut = 0.) //reject clusters with dist to nearest bad channel < cut (in cm)
{
  //Add a task AliAnalysisTaskTaggedPhotons to the analysis train
  //Author: Dmitri Peresunko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTagging", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTagging", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskTaggedPhotons* task = new AliAnalysisTaskTaggedPhotons(Form("%sTask", name));

  task->SelectCollisionCandidates(offlineTriggerMask);

 
  task->SetDistanceToBad() ;
  task->SetTimeCut(25.e-9) ;

 
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), THashList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
