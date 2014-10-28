//AliAnalysisTaskV0sInJets* AddTaskV0sInJets(TString jetBranchName = "", TString outputFile = "output.root", Bool_t bIsMC, TString label = "", Bool_t bTreeOutput = 0, TString outputFileTree = "TreeV0.root")
AliAnalysisTaskV0sInJets* AddTaskV0sInJets(TString jetBranchName = "", TString outputFile = "output.root", Bool_t bIsMC, TString label = "")
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
  {
    Error("AddTaskV0sInJets", "No analysis manager found.");
    return 0;
  }

  TString taskName = "V0";
  TString containerName = "V0histo";
  if(jetBranchName.Length())
  {
    taskName += Form("_%s", jetBranchName.Data());
    containerName += Form("_%s", jetBranchName.Data());
  }
  if(label.Length())
  {
    taskName += Form("_%s", label.Data());
    containerName += Form("_%s", label.Data());
  }
  AliAnalysisTaskV0sInJets* mytask = new AliAnalysisTaskV0sInJets(taskName.Data());

  // Configure task
  mytask->SetJetBranchName(jetBranchName.Data());
  mytask->SetMCAnalysis(bIsMC);
//  mytask->SetTreeOutput(bTreeOutput);

  // Add task
  mgr->AddTask(mytask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* coutput0 = mgr->GetCommonOutputContainer();
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Std"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "QA"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput3 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput4 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "MC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
//  if (bTreeOutput)
//    AliAnalysisDataContainer* coutput5 = mgr->CreateContainer(Form("%s_%s",containerName.Data(),"Tree"), TTree::Class(),AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFileTree.Data(),taskName.Data()));

  // Connect input/output
  mgr->ConnectInput(mytask, 0, cinput0);
  mgr->ConnectOutput(mytask, 0, coutput0);  // No need to connect to a common AOD output container if the task does not fill AOD info.
  mgr->ConnectOutput(mytask, 1, coutput1);
  mgr->ConnectOutput(mytask, 2, coutput2);
  mgr->ConnectOutput(mytask, 3, coutput3);
  mgr->ConnectOutput(mytask, 4, coutput4);
//  if (bTreeOutput)
//    mgr->ConnectOutput(mytask, 5, coutput5);

  return mytask;
}
