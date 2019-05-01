class AliAnalysisDataContainer;
class TNamed;
AliAnalysisTaskGFWFlow* AddTaskGFWFlow(TString name = "name", Bool_t ProduceWeights=kTRUE, Bool_t IsMC=kFALSE, Bool_t AddQA=kFALSE, TString weightpath="",TString subfx="")
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

  //Connect weights to a container
  if(!ProduceWeights) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("InputWeights")) {
      printf("InputWeights not loaded yet, loading now!\n");
      if(weightpath.EqualTo("")) { printf("Weight path for containers not set!\n"); return NULL; };
      if(weightpath.Contains("alien:")) TGrid::Connect("alien:");
      TFile *tf = TFile::Open(weightpath.Data());
      if(!tf) { printf("Could not open weight file %s!\n",weightpath.Data()); return NULL; };
      TList *tl = (TList*)tf->Get("WeightList");
      if(!tl) { printf("Could not wetch WeightList from %s!\n",weightpath.Data()); tf->ls(); return NULL; };
      AliAnalysisDataContainer *cInWeights = mgr->CreateContainer(Form("InputWeights"),TList::Class(), AliAnalysisManager::kInputContainer);
      cInWeights->SetData(tl);
      mgr->ConnectInput(task,1,cInWeights);
    } else {
      mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("InputWeights"));
      printf("InputWeights already loaded\n");
    };
  };

  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1;
  if(ProduceWeights) 
    cOutput1 = mgr->CreateContainer("OutputList", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  else cOutput1 = mgr->CreateContainer(Form("OutCont%s",subfx.Data()), AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  if(AddQA) {
    AliAnalysisDataContainer *qaOutput = mgr->CreateContainer(Form("OutContQA%s",subfx.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
    mgr->ConnectOutput(task,2,qaOutput);
  };


 
  return task;
}
