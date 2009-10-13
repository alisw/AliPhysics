AliAnalysisTaskSECharmFraction* AddTaskCharmFraction(const char* fileout="d0D0.root",Int_t switchMC[5])
{  
  //
  // Configuration macro for the task to analyze the fraction of prompt charm
  // using the D0 impact parameter
  // andrea.rossi@ts.infn.it
  //
  //==========================================================================

  //######## !!! THE SWITCH FOR MC ANALYSIS IS NOT IMPLEMENTED YET!!! ##########à
  switchMC[0]=1;
  switchMC[1]=1;
  switchMC[2]=1;
  switchMC[3]=1;
  switchMC[4]=1;
  Int_t last=0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCharmFraction", "No analysis manager to connect to.");
    return NULL;
  }   

  TString str=fileout,containername;
  str.ReplaceAll(".root","");
  str.Prepend("_");

  AliAnalysisTaskSECharmFraction *hfTask;
 
  hfTask = new AliAnalysisTaskSECharmFraction("AliAnalysisTaskSECharmFraction");
    

  /*  ############### HERE THE POSSIBILITY TO SWITCH ON/OFF THE TLISTS AND MC SELECTION WILL BE SET #########à

  hfTask->SetUseCuts(setD0usecuts);
  hfTask->SetCheckMC(setcheckMC);
  hfTask->SetCheckMC_D0(setcheckMC_D0);
  hfTask->SetCheckMC_2prongs(setcheckMC_2prongs);
  hfTask->SetCheckMC_prompt(setcheckMC_prompt);
  hfTask->SetCheckMC_fromB(setcheckMC_fromB);
  hfTask->SetCheckMC_fromDstar(setSkipD0star);
  hfTask->SetStudyPureBackground(setStudyPureBack);*/
  //  hfTask->SetSideBands(0);
  //  hfTask->SetDebugLevel(2);
  mgr->AddTask(hfTask);
 
 
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput =   mgr->GetCommonInputContainer();
  //mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,cinput);
  

  //Now container for general properties histograms
  containername="coutputNentries";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  
  mgr->ConnectOutput(hfTask,1,coutputNentries);

  containername="coutputSignalType";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  
  mgr->ConnectOutput(hfTask,2,coutputSignalType);


  containername="coutputSignalType_LsCuts";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType_LsCuts = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  
  mgr->ConnectOutput(hfTask,3,coutputSignalType_LsCuts);


 containername="coutputSignalType_TghCuts";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType_TghCuts = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  
  mgr->ConnectOutput(hfTask,4,coutputSignalType_TghCuts);

  // Now container for TLists 
  last=5;
  //##########  NO CUTS TLISTS CONTAINER ##############à
  containername="coutput_nc_sign";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_nc_sign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_nc_sign);
  last++;


  containername="coutput_nc_back";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_nc_back = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_nc_back);
  last++;

  containername="coutput_nc_fromB";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_nc_fromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_nc_fromB);
  last++;


  containername="coutput_nc_fromDstar";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_nc_fromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_nc_fromDstar);
  last++;


  containername="coutput_nc_other";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_nc_other = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_nc_other);
  last++;


  //######### LOOSE CUTS TLISTS CONTAINER #############
  containername="coutput_ls_sign";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_ls_sign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_ls_sign);
  last++;


  containername="coutput_ls_back";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_ls_back = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_ls_back);
  last++;

  containername="coutput_ls_fromB";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_ls_fromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_ls_fromB);
  last++;


  containername="coutput_ls_fromDstar";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_ls_fromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_ls_fromDstar);
  last++;


  containername="coutput_ls_other";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_ls_other = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_ls_other);
  last++;



  //######### TIGHT CUTS TLISTS CONTAINER #############
    containername="coutput_tgh_sign";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_tgh_sign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_tgh_sign);
  last++;


  containername="coutput_tgh_back";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_tgh_back = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_tgh_back);
  last++;

  containername="coutput_tgh_fromB";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_tgh_fromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_tgh_fromB);
  last++;


  containername="coutput_tgh_fromDstar";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_tgh_fromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_tgh_fromDstar);
  last++;


  containername="coutput_tgh_other";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutput_tgh_other = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,last,coutput_tgh_other);
  


  return hfTask;
}
