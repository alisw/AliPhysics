AliAnalysisTaskSEJPSItoEle *AddTaskJPSItoEle(char* fileout = "AliAOD.Jpsitoele.root") 
{
  //***********************************************************************************************
  // Test macro for the AliAnalysisTaskSEJPSItoEle 
  // Starting from AliAOD.root + AliAOD.VertexingHF.root with HF + like sign candidates,
  // it produces a specific AliAODjpsi.root file with a replica of J/psi->e+e- candidates only 
  // (and references to the corresponding decay tracks) + several histograms for
  // both unlike sign and like sign candidates.
  // 
  // C.Di Giglio, carmelo.digiglio@ba.infn.it
  //***********************************************************************************************

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBtoJPSItoEle", "No analysis manager to connect to.");
    return NULL;
  }   

  AliAODHandler* aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName(fileout);

  mgr->SetOutputEventHandler(aodHandler);
  mgr-> SetDebugLevel(10);

  AliAnalysisTaskSEJPSItoEle *hJPSItoEleTask = new AliAnalysisTaskSEJPSItoEle("AOD_JPSItoEle_filter");
  hJPSItoEleTask->SetDebugLevel(2);

  Double_t ptCuts[2] = {0.,500.}; // the cut is on the J/psi pT (--> change this)  
  hJPSItoEleTask->SetPtCuts(ptCuts);
  hJPSItoEleTask->SetAODMCInfo(kFALSE); // only for sim 

  mgr->AddTask(hJPSItoEleTask);

  mgr->ConnectInput(hJPSItoEleTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("tree",TTree::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "default");
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histos",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,       
							   "JPSItoEleHistos.root");
  mgr->ConnectOutput(hJPSItoEleTask,0,coutput0);
  mgr->ConnectOutput(hJPSItoEleTask,1,coutput1);
  return hJPSItoEleTask;
}
