AliAnalysisTaskStrVsMult_RStudies *AddTaskStrVsMult_RStudies(bool K0s = true, bool Lam = true, bool Xi = true, bool Om = true, TString suffix = "")
{
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { 
    ::Error("AddTaskStrangenessVsMultiplicity", "No analysis manager to connect to."); 
    return NULL; 
  }
  if (!mgr->GetInputEventHandler()) { 
    ::Error("AddTaskStrangenessVsMultiplicity", "This task requires an input event handler"); 
    return NULL; 
  }

  // Create the task and add it to the manager
  TString tskname = Form("StrVsMult_Task_%s", suffix.Data());
  AliAnalysisTaskStrVsMult_RStudies *mytask = new AliAnalysisTaskStrVsMult_RStudies(tskname);
  mgr->AddTask(mytask);
  mytask->SetParticleAnalysisStatus(K0s, Lam, Xi, Om);

  // output file name
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGLF_StrVsMult";
  printf("Set OutputFileName : \n %s\n", outputFileName.Data());

  //create and link only used containers
  AliAnalysisDataContainer *coutput[8];
  TString clabels[8] = {"eve", "K0S", "Lam", "ALam", "Xim", "Xip", "Omm", "Omp"};

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer(Form("chists_%s_%s", clabels[0].Data(), suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);mgr->ConnectOutput(mytask, 1, coutput[0]);
  
  int cnumber = 1;
  for (int icont=1; icont<8; icont++) {
    if (mytask->GetParticleAnalysisStatus(icont-1)) {
      coutput[cnumber] = mgr->CreateContainer(Form("chists_%s_%s", clabels[icont].Data(), suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
      mgr->ConnectOutput(mytask, cnumber+1, coutput[cnumber]);
      cnumber++;
    }
  }

  return mytask;

}
