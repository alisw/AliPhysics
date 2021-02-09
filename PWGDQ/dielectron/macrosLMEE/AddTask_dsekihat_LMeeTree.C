AliAnalysisTask *AddTask_dsekihat_LMeeTree(
  const UInt_t trigger = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral,
  const Float_t minPt     = 0.15,
  const Float_t maxEta    = 0.9,
  const Float_t minNsigma = -5.0,
  const Float_t maxNsigma = +5.0
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTask_dsekihat_ReducedTree", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_dsekihat_ReducedTree", "This task requires an input event handler");
    return NULL;
  }

  TString triggername = "NULL";
  if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
  else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
  else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";

  AliAnalysisTaskReducedTreeDS *task = new AliAnalysisTaskReducedTreeDS(Form("ReducedTreeDS_%s",triggername.Data()));
  task->SelectCollisionCandidates(trigger); 
  task->SetMinPtCut(minPt);
  task->SetMaxEtaCut(maxEta);
  task->SetMinTPCNsigmaEleCut(minNsigma);
  task->SetMaxTPCNsigmaEleCut(maxNsigma);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  TString outputFile = AliAnalysisManager::GetCommonFileName();
  //TString listname = "list_BasicInfo";
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname, THashList::Class(), AliAnalysisManager::kOutputContainer,outputFile);
  //mgr->ConnectOutput(task, 1, coutput1);

  const TString treename = "EventTree";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(treename, TTree::Class(), AliAnalysisManager::kOutputContainer,outputFile);
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

