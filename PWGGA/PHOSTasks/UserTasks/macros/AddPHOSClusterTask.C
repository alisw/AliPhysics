#include "/home/alexej/alice/ali-master/AliPhysics/PWGGA/PHOSTasks/UserTasks/AliAnalysisPHOSClusterTask.h"

AliAnalysisPHOSClusterTask* AddPHOSClusterTask (Float_t fZVertex, TString tskname = "AliPHOSCluster") {
cout << "Hello World from the AddTask"<< endl;
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPHOSCluster", "No analysis manager found.");
    return 0x0;
  }

  AliAnalysisPHOSClusterTask* task = new AliAnalysisPHOSClusterTask(tskname.Data());
//  task->SetZVertexCut(fZVertex);
  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("AliPHOSCluster", TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            "AliPHOSCluster.root");



  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("%s_summary", tskname.Data()), TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "AnalysisResults.root");


  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);
  return task;
}
