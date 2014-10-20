AliAnalysisTaskFemtoESE *AddTaskFemtoESE(Double_t qmin = 0, Double_t qmax = 100, Int_t EPdet = 1, Int_t qdet = 0)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFemtoESEtest", "No analysis manager to connect to.");
    return NULL;
  }  

  TString taskname = Form("FemtoESETask_qperc%.2i_%.2i_EP%i_q%i",(Int_t)qmin,(Int_t)qmax,EPdet,qdet);

  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskFemtoESE* ana = new  AliAnalysisTaskFemtoESE(taskname);

  ana->SetQPercCuts(qmin,qmax);
  ana->SetEPDetector(EPdet); // detector used for event plane (0-V0A, 1-V0C)
  ana->SetQPercDetector(qdet); // detector used for q-vector (0-V0A, 1-V0C)

  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  
  TString contname(taskname);
  contname += "_output";

  TString outputFileName = AliAnalysisManager::GetCommonFileName();  

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputA = mgr->CreateContainer(contname, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HelperPIDOutput_%s",contname.Data()), AliHelperPID::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("EventCutsOutput_%s",contname.Data()), AliSpectraAODEventCuts::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("TrackCutsOutput_%s",contname.Data()), AliSpectraAODTrackCuts::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  
  //connect containers
  mgr->ConnectInput(ana, 0, cinput);
  mgr->ConnectOutput(ana, 1, coutputA);
  mgr->ConnectOutput(ana, 2, coutput1);
  mgr->ConnectOutput(ana, 3, coutput2);
  mgr->ConnectOutput(ana, 4, coutput3);
  mgr->AddTask(ana);
  return ana;

}
