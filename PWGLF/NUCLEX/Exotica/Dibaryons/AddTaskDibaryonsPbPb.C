AliAnalysisTask *AddTaskDibaryonsPbPb( Int_t  collidingSystem    = 2,
                                       UInt_t triggerclass       = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral,
                                       Bool_t pileupCut          = kTRUE,
                                       Bool_t pairCleaning       = kFALSE,
                                       Bool_t eventMixing        = kFALSE,
                                       Double_t nsigProton       = 3.0,
                                       Double_t nsigV0Daughter   = 3.0,
                                       Double_t nsigCascDaughter = 3.0 )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskDibaryonsPbPb", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskDibaryonsPbPb", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  AliAnalysisTaskDibaryons *task = new AliAnalysisTaskDibaryons("Dibaryons_PbPb");
  task->SelectCollisionCandidates   (triggerclass); 
  task->SetAnalysisType             (type);
  task->SetCollidingSystem          (collidingSystem);
  task->SetSelectedTriggerClass     (triggerclass);
  task->SetFilterBit                (32);
  task->SetPileupCut                (pileupCut);
  task->SetPairCleaning             (pairCleaning);
  task->SetEventMixing              (eventMixing);
  task->SetNsigmaProton             (nsigProton);
  task->SetNsigmaV0Daughter         (nsigV0Daughter);
  task->SetNsigmaCascDaughter       (nsigCascDaughter);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  TString outputFile = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("list_dibaryon_PbPb", THashList::Class(), AliAnalysisManager::kOutputContainer,outputFile);
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

