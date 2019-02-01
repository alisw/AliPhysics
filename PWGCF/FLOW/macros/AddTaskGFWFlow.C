class AliAnalysisDataContainer;
class TNamed;
AliAnalysisTaskGFWFlow* AddTaskGFWFlow(TString name = "name", Bool_t ProduceWeights=kTRUE, Bool_t IsMC=kFALSE, Bool_t AddQA=kFALSE, TString subfx="")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  if(IsMC) {
    if(!mgr->GetMCtruthEventHandler()) {
      Error("AddTaskGFWFlow","Could not get MC truth handler");
      return NULL;
    };
    AliMCEventHandler *handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    handler->SetReadTR(kTRUE);
  };
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskGFWFlow* task = new AliAnalysisTaskGFWFlow(Form("%s%s",name.Data(),subfx.Data()), ProduceWeights, IsMC, AddQA);
  if(!task)
    return 0x0;
  //My settings:
  mgr->AddTask(task); // add your task to the manager

  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1;
  if(ProduceWeights) 
    cOutput1 = mgr->CreateContainer("OutputList", TList::Class(), AliAnalysisManager::kOutputContainer, "OutputResults.root");
  else cOutput1 = mgr->CreateContainer(Form("OutCont%s",subfx.Data()), AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, Form("OutputResults%s.root",subfx.Data()));
  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  if(AddQA) {
    AliAnalysisDataContainer *qaOutput = mgr->CreateContainer(Form("OutContQA%s",subfx.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("QAResults%s.root",subfx.Data()));
    mgr->ConnectOutput(task,2,qaOutput);
  };
  return task;
}
