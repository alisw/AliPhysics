AliAnalysisTaskSEVertexingHF *AddTaskVertexingHF(Int_t collisionSystem,TString localdir="",TString configfilename="",Int_t runnumber=-1, TString strPeriod="",const char* fname="AliAOD.VertexingHF.root") {
  //
  // Creates a task for heavy flavour vertexing and adds it to the analysis manager.
  // andrea.dainese@lnl.infn.it
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskVertexingHF", "No analysis manager to connect to.");
    return NULL;
  }   
   
  // This task requires an ESD or AOD input handler and an AOD output handler.
  // Check this using the analysis manager.
  //===============================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskVertexingHF", "HF vertexing task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }   
  // Check if AOD output handler exist.
  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskVertexingHF", "HF vertexing task needs the manager to have an AOD output handler.");
    return NULL;
  }

  TString ConfigMode = "";
     if(gSystem->Getenv("CONFIG_MODE"))ConfigMode = gSystem->Getenv("CONFIG_MODE");
  
  // Copy the needed Config file in the current directory 
  if(configfilename.IsNull()){
    TString     configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF.C";// default value: file for pp collisions
    if(collisionSystem==1){
      configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LooseIP.C";
      if(runnumber>=295424)configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt.C";
      if(ConfigMode.Contains("Run3")) configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/upgrade/ConfigVertexingHF_Pb_Upgrade2018.C";
      if(ConfigMode.Contains("ITS3")) configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/upgrade/ConfigVertexingHF_Pb_ITS3.C";
    }
    else if(collisionSystem!=0){
      ::Error("AddTaskVertexingHF","Value of collision system not valid");
    }
    TFile::Cp(gSystem->ExpandPathName(configPWG3d2h.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));
  }
  else{
    TFile::Cp(gSystem->ExpandPathName(configfilename.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));
  }

  // Create the task, add it to the manager and configure it.
  //===========================================================================
  AliAnalysisTaskSEVertexingHF *hfTask = new AliAnalysisTaskSEVertexingHF("vertexing HF");
  hfTask->SetDeltaAODFileName(fname);
  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *coutputListOfCuts = mgr->CreateContainer("ListOfCuts",TList::Class(),AliAnalysisManager::kOutputContainer,hfTask->GetDeltaAODFileName()); //cuts

  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,0,mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(hfTask,1,coutputListOfCuts);

  return hfTask;
}
