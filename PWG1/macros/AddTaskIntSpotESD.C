AliAnalysisTaskIPInfo* AddTaskIntSpotESD() 
{
  //
  // Task to extract the Int.Spot position and sigma as well as the vertex and track 
  // DCA resolutions. Performs estimates both with TPC and TPC+ITS tracks. 
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskIntSpotESD", "No analysis manager to connect to.");
    return 0;
  }   
  //

  // Create the task
  AliAnalysisTaskIPInfo *taskIP = new AliAnalysisTaskIPInfo("IPInfo");
  taskIP->SetOptions(AliAnalysisTaskIPInfo::kITSTPC, kFALSE, 1e-4, 2, 12,1000,
		     -4e-2,8e-2,  10,2,32,  14,0.2,3., kFALSE);
  taskIP->SetIPCenIni(AliAnalysisTaskIPInfo::kITSTPC, -0.0764,0.2481,0);
  //
  //
  taskIP->SetOptions(AliAnalysisTaskIPInfo::kTPC, kFALSE, 1e-4,	2, 12,1000,
		     -4e-2,8e-2,  10,2,32,  14,0.2,3., kFALSE);
  taskIP->SetIPCenIni(AliAnalysisTaskIPInfo::kTPC, -0.0764,0.2481,0);
  //
  //
  taskIP->SetOptions(AliAnalysisTaskIPInfo::kSPD, kFALSE, 1e-4,	2, 12,1000,
		     -4e-2,8e-2,  10,2,32,  14,0.2,3., kTRUE);
  taskIP->SetIPCenIni(AliAnalysisTaskIPInfo::kSPD, -0.0764,0.2481,0);
  //
  mgr->AddTask(taskIP);

  // Create containers for input/output
  AliAnalysisDataContainer *cInputIPesd = mgr->CreateContainer
    ("cInputIPesd",TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *cOutputIPesd = mgr->CreateContainer
    ("cOutputIPesd",TList::Class(),AliAnalysisManager::kOutputContainer,"IPInfo.root");

  // Attach input
  mgr->ConnectInput(taskIP,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskIP,0,cOutputIPesd);
  
  return taskIP;
}
