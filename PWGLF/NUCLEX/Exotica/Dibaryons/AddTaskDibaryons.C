AliAnalysisTask *AddTaskDibaryons( Int_t   collidingSystem                      = 0,
                                   AliVEvent::EOfflineTriggerTypes triggerclass = AliVEvent::kINT7,
                                   Bool_t  pileupCut                            = kTRUE,
                                   Bool_t  eventMixing                          = kTRUE )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskDibaryons", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskDibaryons", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  TString triggername = "NULL";
  if(triggerclass == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
  else if(triggerclass == (UInt_t)AliVEvent::kHighMultV0)  triggername = "kHighMultV0";
  else if(triggerclass == (UInt_t)AliVEvent::kHighMultSPD) triggername = "kHighMultSPD";

  AliAnalysisTaskDibaryons *task = new AliAnalysisTaskDibaryons(Form("Dibaryons_%s",triggername.Data()));
  task->SelectCollisionCandidates   (triggerclass); 
  task->SetAnalysisType             (type);
  task->SetCollidingSystem          (collidingSystem);
  task->SetSelectedTriggerClass     (triggerclass);
  task->SetFilterBit                (128);
  task->SetPileupCut                (pileupCut);
  task->SetEventMixing              (eventMixing);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  TString outputFile = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("list_dibaryon_%s",triggername.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer,outputFile);
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

