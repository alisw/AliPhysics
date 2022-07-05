
AliAnalysisTaskElectronEfficiencyV2* AddTask_sscheid_efficiencyV2(TString name = "name",
                                                                Bool_t getFromAlien = kFALSE,
                                                                TString configFile="Config_sscheid_efficiencyV2.C",
                                                                const std::string resolutionFilename ="",
                                                                Int_t wagonnr = 0
                                                                ) {

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  // #########################################################
  // #########################################################
  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName = "AnalysisResults.root"; // create a subfolder in the file

  // #########################################################
  // #########################################################
  // Loading individual config file either local or from Alien

  // TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //Load updated macros from private ALIEN path
  if (!gSystem->AccessPathName(configFile))
  {
    printf("file already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
  else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hscheid/PWGDQ/dielectron/macrosLMEE/%s file:./",configFile.Data()))) )
  {
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+configFile);
  // #########################################################
  // #########################################################
  // Creating an instance of the task

  //Get instance of task from config, get rif of cutlib
  TMacro conf_die(gSystem->ExpandPathName(configFilePath.Data())); //ROOT6
  AliAnalysisTaskElectronEfficiencyV2 *task = reinterpret_cast<AliAnalysisTaskElectronEfficiencyV2 *>(conf_die.Exec(Form("\"%s\",%d",name.Data(),wagonnr)));

  // add post calibration files to task
  task->SetResolutionFile(resolutionFilename,"/alice/cern.ch/user/h/hscheid/supportFiles/" + resolutionFilename);



  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d", wagonnr), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
