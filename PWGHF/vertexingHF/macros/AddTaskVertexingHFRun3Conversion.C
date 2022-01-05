AliAnalysisTaskSEVertexingHF *AddTaskVertexingHFRun3Conversion(TString configfilename=""){
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

  TString localdir=".";
  Int_t runnumber=-1;
  if(gSystem->Getenv("ALIEN_JDL_LPMRUNNUMBER")) runnumber = atoi(gSystem->Getenv("ALIEN_JDL_LPMRUNNUMBER"));
  TString collSyst=gSystem->Getenv("ALIEN_JDL_LPMINTERACTIONTYPE");
  TString prodType = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
  
  // Copy the needed Config file in the current directory 
  if(configfilename.IsNull()){
    TString     configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF.C";
    if(collSyst=="PbPb" || collSyst=="XeXe"){
      configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LooseIP.C";
      if(runnumber>=295424)configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt.C";
      if (prodType=="MC"){
	TString prodtag=gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTAG");
	if(prodtag=="LHC20g2a_2" || prodtag=="LHC20g11a2" || prodtag=="LHC20j5a1" || prodtag=="LHC20j5a2" || prodtag=="LHC20k3a" || prodtag=="LHC20l5a") configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt_Central.C";
	else if(prodtag=="LHC20g2b_2" || prodtag=="LHC20g11b2" || prodtag=="LHC20j5b1" || prodtag=="LHC20j5b2" || prodtag=="LHC20k3b" || prodtag=="LHC20l5b") configPWG3d2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LcMinpt1_DsMinPt15_2018opt_SemiCentral.C";
      }
    }
    Printf("HF config file that will be used (and copied to ConfigVertexingHF.C in the local directory) is: %s", configPWG3d2h.Data());
    TFile::Cp(gSystem->ExpandPathName(configPWG3d2h.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));
  }else{
    Printf("HF config file that will be used (and copied to ConfigVertexingHF.C in the local directory) is: %s", configfilename.Data());
    TFile::Cp(gSystem->ExpandPathName(configfilename.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));
  }

  // Create the task, add it to the manager and configure it.
  //===========================================================================
  AliAnalysisTaskSEVertexingHF *hfTask = new AliAnalysisTaskSEVertexingHF("vertexing HF");
  hfTask->SetDeltaAODFileName("AliAOD.VertexingHF.root");
  mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *coutputListOfCuts = mgr->CreateContainer("ListOfCuts",TList::Class(),AliAnalysisManager::kOutputContainer,hfTask->GetDeltaAODFileName()); //cuts

  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,0,mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(hfTask,1,coutputListOfCuts);

  return hfTask;
}
