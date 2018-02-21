AliAnalysisTaskTOFMC *AddTaskTOFMC(const char *PeriodName=NULL, Int_t nTPC_CR=70, Int_t Chi2_TPCcluser=4, Int_t DCAz=2)
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskTOFMC", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskTOFMC", "This task requires an input event handler");
      return NULL;
   } 
	
	TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddTaskTOFMC", "input data type is not \"ESD\"");
    return NULL;
  }

  /* check MC truth event handler */
    if (!mgr->GetMCtruthEventHandler()) {
      Error("AddTaskTOFMC", "cannot get MC truth event handler");
      return NULL;
    }
    AliMCEventHandler * handler = (AliMCEventHandler *) mgr->GetMCtruthEventHandler();
    handler->SetReadTR(kTRUE);


	AliMultSelectionTask *mult = (AliMultSelectionTask*)mgr->GetTask("taskMultSelection");
      mult->SetUseDefaultMCCalib(kTRUE);
      mult->SetAlternateOADBforEstimators(PeriodName);
 
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create and configure the task
  AliAnalysisTaskTOFMC *taskTOFMC = new AliAnalysisTaskTOFMC(PeriodName, nTPC_CR,Chi2_TPCcluser, DCAz);

	AliESDtrackCuts *fTrackCuts =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
	fTrackCuts->SetMinNCrossedRowsTPC(nTPC_CR);
	fTrackCuts->SetMaxChi2PerClusterTPC(Chi2_TPCcluser);
	fTrackCuts->SetMaxDCAToVertexZ(DCAz);
  taskTOFMC->SetTrackCuts(fTrackCuts);
  taskTOFMC->SetTrackCuts2(fTrackCuts);

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
