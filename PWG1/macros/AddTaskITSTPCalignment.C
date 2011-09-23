AliAnalysisTaskITSTPCalignment *AddTaskITSTPCalignment()
{
  //add the ITS TPC alignemtn task to the manager
  //Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch

  //______________________________________________________________________________
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskITSTPCalignment", "No analysis manager to connect to.");
    return NULL;
  }

  //______________________________________________________________________________
  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskITSTPCalignment", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type!="ESD")
  {
    ::Error("AddTaskITSTPCalignment", "This task only works on ESDs");
    return NULL;
  }

  //___________________________________________________________________________
  // Create the task, add it to manager and configure it.
  AliAnalysisTaskITSTPCalignment *task = new AliAnalysisTaskITSTPCalignment("taskITSTPCalignment");
  TTimeStamp t0(2009,11,1,0,0,0);
  TTimeStamp tend(2012,12,31,0,0,0);
  Int_t slotwidth = 3600;
  task->SetupAlignerArray(t0.GetSec(),tend.GetSec(),slotwidth);
  task->SetFillDebugTree(kFALSE);
  task->SetDoQA(kTRUE);
  task->SetMinPt(0.4);
  task->SetMinNclsITS(4);
  task->SetMinNclsTPC(70);
  task->SetRejectOutliers(kTRUE); //internal KF outlier rejection (kalman update-based)
  task->SetRejectOutliersSigma2Median(kTRUE); //input data outlier removal
  task->SetOutRejSigma(1.); //max distance the kf state is allowed to jump
  task->SetOutRejSigmaOnMerge(10.); //outlier rejection when merging vertically
  task->SetOutRejSigma2Median(2.); //max distance from median for input data
  task->SetUseITSoutGlobalTrack(kFALSE);
  task->SetUseITSoutITSSAtrack(kTRUE);
  
  mgr->AddTask(task);

  //______________________________________________________________________________
  //connect output
  TString outputFilename = "outputITSTPCalignment.root";

  AliAnalysisDataContainer* coutput0 = mgr->CreateContainer("outputTree",
                                       TTree::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputFilename.Data());
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("outputList",
                                       TList::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputFilename.Data());
  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer("outputArrayITSglobal",
                                       AliRelAlignerKalmanArray::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputFilename.Data());
  AliAnalysisDataContainer* coutput3 = mgr->CreateContainer("outputArrayITSSA",
                                       AliRelAlignerKalmanArray::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputFilename.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);

  return task;
}

