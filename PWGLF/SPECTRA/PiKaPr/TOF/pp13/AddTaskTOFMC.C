AliAnalysisTaskTOFMC *AddTaskTOFMC(Int_t nTPC_CR=70, Int_t Chi2_TPCcluser=4, Int_t DCAz=2)
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskStrangenessVsMultiplicity", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskdE_D", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create and configure the task
  AliAnalysisTaskTOFMC *taskTOFMC = new AliAnalysisTaskTOFMC("AliAnalysisTaskTOFMC", nTPC_CR,Chi2_TPCcluser, DCAz);

	AliESDtrackCuts *fesdTrackCuts =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
	fesdTrackCuts->SetMinNCrossedRowsTPC(nTPC_CR);
	fesdTrackCuts->SetMaxChi2PerClusterTPC(Chi2_TPCcluser);
	fesdTrackCuts->SetMaxDCAToVertexZ(DCAz);

  taskTOFMC->SetTrackCuts(fesdTrackCuts);

  mgr->AddTask(taskTOFMC);
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
  outputFileName += ":PWGLF_pp13VsMult";
  TString OutputListname = Form("fOutputList_CR%i_Chi2TPCcluser%i_DCAz%i",nTPC_CR, Chi2_TPCcluser, DCAz);
//  if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   
  Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

  AliAnalysisDataContainer *coutputList = mgr->CreateContainer(OutputListname,
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, //"MC_PiKaPr_TOF_std.root");
							     outputFileName );

  mgr->ConnectInput (taskTOFMC, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskTOFMC, 1, coutputList);
  
  return taskTOFMC;
}
