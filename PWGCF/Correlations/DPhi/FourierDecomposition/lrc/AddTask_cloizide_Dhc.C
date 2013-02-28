AliAnalysisTask *AddTask_cloizide_Dhc() 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_cloizide_Dhc", "No analysis manager found.");
    return 0;
  }

  TString list=gSystem->Getenv("LIST");
  Bool_t isPass1 = list.Contains("pass1");
  Bool_t isPass2 = list.Contains("pass2");

  if (isPass1)
    return 0;

  gROOT->ProcessLine(".L KiddiePoolClasses.cxx+");
  gROOT->ProcessLine(".L AliDhcTask_opt.cxx+");

  AliDhcTask *dhcTask = new AliDhcTask("Task_cloizide_Dhc");
  dhcTask->SelectCollisionCandidates(AliVEvent::kMB);
  dhcTask->SetVerbosity(10);
  mgr->AddTask(dhcTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(dhcTask,0,cinput);
  AliAnalysisDataContainer *co1 = mgr->CreateContainer("Cont_cloizide_DhcAna",
                                                       TList::Class(), 
                                                       AliAnalysisManager::kOutputContainer,
                                                       "cloizide_DhcAna.root");
  mgr->ConnectOutput(dhcTask,1,co1);

  return 0;
}
