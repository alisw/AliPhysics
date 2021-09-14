AliAnalysisTaskBEpp13TeV* AddTaskBEpp13TeV(
	TString	taskName = "beauty",
	bool	isMC = false,
	int		minTPCnCrossedRow = 100,
	int		minTPCnClsPID = 80,
	double	maxTPCchi2 = 4,
	double	minTPCclsRatio = 0.6,
	int		minITSnCls = 3,
	int		itsLayer = 2,
	double	tpcPIDlow = -1,
	double	tpcPIDhigh = 3,
	double	tofPID = 3)
{
  //TString taskName = "beauty";
  //TString itsLayer = citsLayer;
  printf("=======================================\n");
  printf("==========ADD TASK PARAMETERS==========\n");
  printf("=======================================\n");
  printf("Task name: %s \n", taskName.Data());
  printf("MC Flag: %d \n",isMC);
  printf("Min. TPC crossed raw: %d \n", minTPCnCrossedRow);
  printf("Min. TPC clusters for PID: %d \n", minTPCnClsPID);
  printf("Max. TPC chi2/ndf: %.1f \n", maxTPCchi2);
  printf("Min. TPC cls/clsPID: %.1f \n", minTPCclsRatio);
  printf("Min. ITS clusters: %d \n", minITSnCls);
  printf("ITS requirement: %d \n",itsLayer);
  printf("Min. TPC nsigma: %.1f \n", tpcPIDlow);
  printf("Max. TPC nsigma: %.1f \n", tpcPIDhigh);
  printf("TOF nsigma: %.1f \n", tofPID);
  printf("===================================\n");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    return 0x0;
  }
  
  if(!mgr->GetInputEventHandler()){
    return 0x0;
  }
  
  // now we create an instance of your task
  AliAnalysisTaskBEpp13TeV* task = new AliAnalysisTaskBEpp13TeV(taskName.Data());   
  if(!task) return 0x0;
  if(isMC) task->SetMCanalysis();
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  task->SetMinTPCnCrossedRow(minTPCnCrossedRow);
  task->SetMinTPCNclsPID(minTPCnClsPID);
  task->SetMaxTPCchi2(maxTPCchi2);
  task->SetMinTPCclsRatio(minTPCclsRatio);
  task->SetMinITSNcls(minITSnCls);
  task->SetITSlayer(itsLayer);
  task->SetPIDCuts(tpcPIDlow, tpcPIDhigh, tofPID);
  // add your task to the manager
  mgr->AddTask(task);
  
  // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":ElectronID_";      // create a subfolder in the file
  fileName += taskName;

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("ccontainer0_%s",taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  return task;
}
