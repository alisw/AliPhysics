AliAnalysisTaskSEVertexingHFRun3Conversion *AddTaskVertexingHFRun3Conversion(TString configfilename="", Bool_t disableCasc=kFALSE){
  //
  // Creates a task for heavy flavour vertexing for conversion to AO2D
  // Extract parameters from environment variables
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskVertexingHF", "No analysis manager to connect to.");
    return NULL;
  }   
   
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskVertexingHF", "HF vertexing task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }   

  auto getProdData = [](int &runnumber, TString &collSys, TString &prodType, TString &prodTag)
  {
    // This works in the lego trains where the file env.sh is sourced
    runnumber = -1;
    if (gSystem->Getenv("ALIEN_JDL_LPMRUNNUMBER")) {
      runnumber = atoi(gSystem->Getenv("ALIEN_JDL_LPMRUNNUMBER"));
      collSys   = gSystem->Getenv("ALIEN_JDL_LPMINTERACTIONTYPE");
      prodType  = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
      prodTag   = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTAG");
      return 0;
    }
    if (!gSystem->Getenv("CHILD_DATASETS")) {
      ::Error("AddTaskVertexingHF","Environment variable CHILD_DATASETS not found. Cannot configure task.");
      return 1;
    }
    int num_children = atoi(gSystem->Getenv("CHILD_DATASETS"));
    char varname[50];
    for (int i = 1; i <= num_children; ++i) {
      // Check if child is active
      sprintf(varname, "RUNNO_child_%d", i);
      if (!gSystem->Getenv(varname) || atoi(gSystem->Getenv(varname)) < 0) continue;
      sprintf(varname, "ALIEN_JDL_child_%d_LPMRUNNUMBER", i);
      runnumber = atoi(gSystem->Getenv(varname));
      sprintf(varname, "ALIEN_JDL_child_%d_LPMINTERACTIONTYPE", i);
      collSys   = gSystem->Getenv(varname);
      sprintf(varname, "ALIEN_JDL_child_%d_LPMPRODUCTIONTYPE", i);
      prodType  = gSystem->Getenv(varname);
      sprintf(varname, "ALIEN_JDL_child_%d_LPMPRODUCTIONTAG", i);
      prodTag  = gSystem->Getenv(varname);
      return 0;
    }
    // No child is active -> error
    ::Error("AddTaskVertexingHF","No child of this train seems to be active. Cannot configure task.");
    return 1;
  };

  TString localdir=".";
  Int_t runnumber=-1;
  TString collSyst="";
  TString prodType="";
  TString prodTag="";
  Int_t readEnv=getProdData(runnumber, collSyst, prodType, prodTag);
  if(readEnv!=0) return NULL;
  
  ::Info("AddTaskVertexingHF","Found: runno=%d collSys=%s prodType=%s prodTag=%s\n", runnumber, collSyst.Data(), prodType.Data(),prodTag.Data());
  // Copy the needed Config file in the current directory 
  if(configfilename.IsNull()){
    if(collSyst.IsNull()){
      ::Error("AddTaskVertexingHF","LPMInteractionType not available and custom config not passed in the arguments");
      return NULL;
    }
    TString configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_pp_OfflineV0.C";
    if(collSyst=="PbPb" || collSyst=="XeXe"){
      configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LooseIP.C";
      if(runnumber<0){
        ::Error("AddTaskVertexingHF","LPMRunNumber (needed for Pb-Pb or Xe-Xe) not available and custom config not passed in the arguments");
        return NULL;
      }
      if(prodType.IsNull()){
        ::Error("AddTaskVertexingHF","LPMProductionType (needed for Pb-Pb or Xe-Xe) not available and custom config not passed in the arguments");
        return NULL;
      }
      if(runnumber>=295424)configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt.C";
      if (prodType=="MC"){
        if(prodTag.IsNull()){
          ::Error("AddTaskVertexingHF","LPMProductionTag (needed for MCs of Pb-Pb or Xe-Xe) not available and custom config not passed in the arguments");
          return NULL;
        }
        if(prodTag=="LHC20g2a_2" || prodTag=="LHC20g11a2" || prodTag=="LHC20j5a1" || prodTag=="LHC20j5a2" || prodTag=="LHC20k3a" || prodTag=="LHC20l5a") configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt_Central.C";
        else if(prodTag=="LHC20g2b_2" || prodTag=="LHC20g11b2" || prodTag=="LHC20j5b1" || prodTag=="LHC20j5b2" || prodTag=="LHC20k3b" || prodTag=="LHC20l5b") configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt_SemiCentral.C";
      }
    }
    Printf("HF config file that will be used (and copied to ConfigVertexingHF.C in the local directory) is: %s", configPWG3d2h.Data());
    TFile::Cp(gSystem->ExpandPathName(configPWG3d2h.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));
  }else{
    Printf("HF config file passed as argument that will be used (and copied to ConfigVertexingHF.C in the local directory) is: %s", configfilename.Data());
    TFile::Cp(gSystem->ExpandPathName(configfilename.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));
  }

  // Create the task, add it to the manager and configure it.
  //===========================================================================
  AliAnalysisTaskSEVertexingHFRun3Conversion *hfTask = new AliAnalysisTaskSEVertexingHFRun3Conversion("vertexing HF");
  if(disableCasc) hfTask->ForceDisableCascades();
  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutputListOfCuts = mgr->CreateContainer("ListOfCuts",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); //cuts
  
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,1,coutputListOfCuts);

  return hfTask;
}
