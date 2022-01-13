AliAnalysisTaskSEVertexingHFRun3Conversion *AddTaskVertexingHFRun3Conversion(TString configfilename="", Bool_t resetAtEachEv=kTRUE){
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
  AliAnalysisTaskSEVertexingHFRun3Conversion *hfTask = new AliAnalysisTaskSEVertexingHFRun3Conversion("vertexing HF");
  hfTask->ResetTreeAtEachEvent(resetAtEachEv);
  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutputListOfCuts = mgr->CreateContainer("ListOfCuts",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); //cuts
  
  AliAnalysisDataContainer *coutputD0 = mgr->CreateContainer("D0CandidateTree",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputfile.Data());
  coutputD0->SetSpecialOutput();
  
  AliAnalysisDataContainer *coutput3p = mgr->CreateContainer("Charm3pCandidateTree",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputfile.Data());
  
  coutput3p->SetSpecialOutput();
  
  AliAnalysisDataContainer *coutputDst = mgr->CreateContainer("DstarCandidateTree",
							      TTree::Class(),
							      AliAnalysisManager::kOutputContainer,
							      outputfile.Data());
  coutputDst->SetSpecialOutput();

  AliAnalysisDataContainer *coutputCasc = mgr->CreateContainer("LcV0bachCandidateTree",
							       TTree::Class(),
							       AliAnalysisManager::kOutputContainer,
							       outputfile.Data());
  coutputCasc->SetSpecialOutput();


  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,1,coutputListOfCuts);
  mgr->ConnectOutput(hfTask,2,coutputD0);
  mgr->ConnectOutput(hfTask,3,coutput3p);
  mgr->ConnectOutput(hfTask,4,coutputDst);
  mgr->ConnectOutput(hfTask,5,coutputCasc);

  return hfTask;
}
