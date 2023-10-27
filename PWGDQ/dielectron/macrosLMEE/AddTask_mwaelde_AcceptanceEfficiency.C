AliAnalysisTask* AddTask_mwaelde_AcceptanceEfficiency(TString name = "name", Bool_t getFromAlien  = kFALSE, TString cFileName ="Config_mwaelde_AccEff.C", Bool_t sysUncOutput  = kFALSE)
{
  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  if (!gSystem->AccessPathName(cFileName)) {
    printf("Configfile already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
  else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/mwalde/PWGDQ/dielectron/macrosLMEE/%s file:./",cFileName.Data()))) ){
    printf("Copy configfile from alien\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config_mwaelde_AccEff")) {
    printf("Load macro now\n");
    gROOT->LoadMacro(configFilePath.Data());
  }

  // get the input event handler, again via a static method.
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MyTask";      // create a subfolder in the file
  AliAnalysisTaskOmegaDielectron_AccEff* task = new AliAnalysisTaskOmegaDielectron_AccEff(name.Data());
  if(!task) return 0x0;
  // task->SelectCollisionCandidates(AliVEvent::kINT7);// kAnyINT
  task->SetBoolsysUncOutput(sysUncOutput);


  // Number of cuts         // "PWGDQ/dielectron/macrosLMEE/AddTask_rbailhac_ElectronEfficiencyV2_PbPb.C"
  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  // std::cout << "  fFilter_TrackCuts has the SIZE : " << nDie << "  !!  " << '\n';

  //add dielectron analysis with different cuts to the task
  for (Int_t i=1; i<=nDie; i++){ //nDie defined in config file
    AliAnalysisFilter *filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("Config_mwaelde_AccEff(%d)",i)));
    task->AddTrackCut(filter);
    AliAnalysisFilter *filter_PIDCuts = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("Config_mwaelde_AccEff_PID(%d)",i)));
    task->AddPIDCut(filter_PIDCuts);
  }

  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("MyOutputContainer_cutVari_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;

}
