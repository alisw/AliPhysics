class AliAnalysisDataContainer;
class TNamed;
AliAnalysisTaskGFWPIDFlow* AddTaskGFWPIDFlow(TString name = "name", Bool_t IsMC=kFALSE, TString stage = "dev", Int_t GFWMode=0,
                          TString NUAWeights="", TString PIDWeights="alien:///alice/cern.ch/user/v/vvislavi/Weights/PIDWeights.root", TString PIDWeights2="")
{
  Int_t StageSwitch = 0;
  if(stage.Contains("weights")) StageSwitch=1;
  if(stage.Contains("meanpt")) StageSwitch=2;
  if(stage.Contains("full")) StageSwitch=3;
  if(stage.Contains("dev")) StageSwitch=4;
  if(stage.Contains("CustomWeights")) StageSwitch=5;
  if(StageSwitch==0) return 0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  if(IsMC) {
    if(!mgr->GetMCtruthEventHandler()) {
      Error("AddAliAnalysisTaskGFWPIDFlow","Could not get MC truth handler");
      return NULL;
    };
    AliMCEventHandler *handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    handler->SetReadTR(kTRUE);
  };
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskGFWPIDFlow* task = new AliAnalysisTaskGFWPIDFlow(name.Data(), IsMC, stage);
  task->SetGFWMode(GFWMode);
  if(!task)
    return 0x0;

  mgr->AddTask(task); // add your task to the manager

  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task


  //Producing weights
  if(StageSwitch==1) {
    AliAnalysisDataContainer *weightCont = mgr->CreateContainer("WeightList",TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,weightCont);
    return task;
  }
  //Load input mean pt
  if(StageSwitch==2) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("Weights")) {
      TGrid::Connect("alien:");
      TFile *tfWeights = TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MergedWeights.root");
      if(!tfWeights) AliFatal("Could not open weights file\n");
      if(tfWeights->IsZombie()) AliFatal("Weight file is a zombie\n");
      TList *fList = (TList*)tfWeights->Get("WeightList");
      if(!fList) { AliFatal("Could not fetch the weight list!\n"); return; };
      AliAnalysisDataContainer *cWeights = mgr->CreateContainer("Weights",TList::Class(), AliAnalysisManager::kInputContainer);
      cWeights->SetData(fList);
      mgr->ConnectInput(task,1,cWeights);
    };
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer("MPTProfileList", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    AliAnalysisDataContainer *multiDist  = mgr->CreateContainer("multiDist",TH1D::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,2,multiDist);
    return task;
  };
  //Full
  if(StageSwitch==3) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("Weights")) {
      TGrid::Connect("alien:");
      TFile *tfWeights = TFile::Open(NUAWeights.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MergedWeights.root"
      TList *fList = (TList*)tfWeights->Get("WeightList");
      AliAnalysisDataContainer *cWeights = mgr->CreateContainer("WeightList",TList::Class(), AliAnalysisManager::kInputContainer);
      cWeights->SetData(fList);
      mgr->ConnectInput(task,1,cWeights);
    };
    if(!AllContainers->FindObject("PIDWeights")) {
      // TGrid::Connect("alien:");
      TFile *tfMPT = TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MeanPts_05_20.root");
      TList *fMPTList = (TList*)tfMPT->Get("MeanPts");
      AliAnalysisDataContainer *cInMPT = mgr->CreateContainer("InputMeanPt",TList::Class(), AliAnalysisManager::kInputContainer);
      cInMPT->SetData(fMPTList);
      mgr->ConnectInput(task,2,cInMPT);
    };
    AliAnalysisDataContainer *cOutputCOV = mgr->CreateContainer("MPTDiff",TProfile::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputCOV);
    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer("FlowCont",AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,2,cOutputFC);
    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer("Covariance",TProfile::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,3,cOutputFC);
    return task;
  };
  if(StageSwitch==4) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("PIDWeights")) {
      TGrid::Connect("alien:");
      TFile *tfWeights = TFile::Open(PIDWeights.Data());
      TList *tlWeights = (TList*)tfWeights->Get("PIDWeights");
      AliAnalysisDataContainer *cInWeights = mgr->CreateContainer("PIDWeights",TList::Class(), AliAnalysisManager::kInputContainer);
      cInWeights->SetData(tlWeights);
      mgr->ConnectInput(task,1,cInWeights);
    } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("PIDWeights"));
    if(!PIDWeights2.IsNull()) {
      if(!AllContainers->FindObject("PIDWeights_ZM")) {
        TGrid::Connect("alien:");
        TFile *tfWeights = TFile::Open(PIDWeights2.Data());
        TList *tlWeights = (TList*)tfWeights->Get("PIDWeights");
        AliAnalysisDataContainer *cInWeights = mgr->CreateContainer("PIDWeights_ZM",TList::Class(), AliAnalysisManager::kInputContainer);
        cInWeights->SetData(tlWeights);
        mgr->ConnectInput(task,2,cInWeights);
      } else mgr->ConnectInput(task,2,(AliAnalysisDataContainer*)AllContainers->FindObject("PIDWeights_ZM"));

    }

    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer(Form("FlowCont_%i",GFWMode),AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputFC);
    return task;
  };
  if(StageSwitch==5) {
    AliAnalysisDataContainer *cOutputWeights  = mgr->CreateContainer("CustomWeights",TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputWeights);
    return task;
  };

  return 0;
}
