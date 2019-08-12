AliAnalysisTaskTOFppSpectra *AddTaskTOFppSpectra(const char *PeriodName=NULL, Int_t nTPC_CR=70, Int_t Chi2_TPCcluser=4, Int_t DCAz=2, Int_t DCAxy=7, Int_t value_Sigma=80, Int_t value_Slope=0)
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskTOFppSpectra", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskTOFppSpectra", "This task requires an input event handler");
      return NULL;
   } 
	
	TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddTaskTOFppSpectra", "input data type is not \"ESD\"");
    return NULL;
  }

	AliMultSelectionTask *mult = (AliMultSelectionTask*)mgr->GetTask("taskMultSelection");
      mult->SetUseDefaultCalib(kTRUE);
      mult->SetAlternateOADBforEstimators(PeriodName);
 
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create and configure the task
  AliAnalysisTaskTOFppSpectra *taskTOFppSpectra = new AliAnalysisTaskTOFppSpectra(PeriodName, nTPC_CR,Chi2_TPCcluser, DCAz, DCAxy, value_Sigma, value_Slope);

	taskTOFppSpectra->SetMinNCrossedRowsTPC(nTPC_CR);
	taskTOFppSpectra->SetMaxChi2PerClusterTPC(Chi2_TPCcluser);
	taskTOFppSpectra->SetMaxDCAToVertexZ(DCAz);

/*	taskTOFppSpectra->SetParamater(2, value_Sigma);
	taskTOFppSpectra->SetParamater(4, value_Slope;
*/
	if (DCAxy==7) taskTOFppSpectra->SetDCAtoVertexXYPtDep("0.0105+0.0350/pt^1.1");
	if (DCAxy==6) taskTOFppSpectra->SetDCAtoVertexXYPtDep("0.0090+0.0300/pt^1.1");
	if (DCAxy==8) taskTOFppSpectra->SetDCAtoVertexXYPtDep("0.0120+0.0400/pt^1.1");

	mgr->AddTask(taskTOFppSpectra);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
	if (value_Slope==0)	Float_t slope=0.0125;
	else if (value_Slope==+1)	Float_t slope=0.01375;
	else if (value_Slope==-1)	Float_t slope=0.01125;

  outputFileName += ":Data_PWGLF_pp13VsMult";
  TString OutputListname = Form("fOutputList_CR%i_Chi2TPCcluster%i_DCAz%i_DCAxy%iSigma_PID_Sigma%i_Slope%f",nTPC_CR, Chi2_TPCcluser, DCAz, DCAxy, value_Sigma, slope);
//  if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   
  Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

  AliAnalysisDataContainer *coutputList = mgr->CreateContainer(OutputListname,
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, //"MC_PiKaPr_TOF_std.root");
							     outputFileName );

  mgr->ConnectInput (taskTOFppSpectra, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskTOFppSpectra, 1, coutputList);
  
  return taskTOFppSpectra;
}
