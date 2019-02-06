AliAnalysisTaskTOFMC *AddTaskTOFMC(const char *PeriodName=NULL, Int_t nTPC_CR=70, Int_t Chi2_TPCcluser=4, Int_t DCAz=2, Int_t DCAxy=7)
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
  AliAnalysisTaskTOFMC *taskTOFMC = new AliAnalysisTaskTOFMC(PeriodName, nTPC_CR,Chi2_TPCcluser, DCAz, DCAxy);

	taskTOFMC->SetMinNCrossedRowsTPC(nTPC_CR);
	taskTOFMC->SetMaxChi2PerClusterTPC(Chi2_TPCcluser);
	taskTOFMC->SetMaxDCAToVertexZ(DCAz);

	if (DCAxy==7) taskTOFMC->SetDCAtoVertexXYPtDep("0.0105+0.0350/pt^1.1");
	if (DCAxy==6) taskTOFMC->SetDCAtoVertexXYPtDep("0.0090+0.0300/pt^1.1");
	if (DCAxy==8) taskTOFMC->SetDCAtoVertexXYPtDep("0.0120+0.0400/pt^1.1");

	mgr->AddTask(taskTOFMC);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
  outputFileName += ":PWGLF_pp13VsMult";
  TString OutputListname = Form("fOutputList_CR%i_Chi2TPCcluster%i_DCAz%i_DCAxy_%iSigma",nTPC_CR, Chi2_TPCcluser, DCAz, DCAxy);
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
