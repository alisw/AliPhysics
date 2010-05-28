AliAnalysisTaskSECharmFraction* AddTaskCharmFraction(TString fileout="d0D0",Int_t switchMC[5],Bool_t readMC=kTRUE)
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

  TString outfile=AliAnalysisManager::GetCommonFileName();
  outfile += ":PWG3_D2H";
  outfile += str.Data();

  AliAnalysisTaskSECharmFraction *hfTask;
 
  hfTask = new AliAnalysisTaskSECharmFraction("AliAnalysisTaskSECharmFraction");
  hfTask->SetReadMC(readMC);    

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
							   outfile.Data());
  
  mgr->ConnectOutput(hfTask,1,coutputNentries);

  containername="coutputSignalType";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  
  mgr->ConnectOutput(hfTask,2,coutputSignalType);


  containername="coutputSignalType_LsCuts";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType_LsCuts = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  
  mgr->ConnectOutput(hfTask,3,coutputSignalType_LsCuts);


 containername="coutputSignalType_TghCuts";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType_TghCuts = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  
  mgr->ConnectOutput(hfTask,4,coutputSignalType_TghCuts);

  // Now container for TLists 
  last=5;
  //##########  NO CUTS TLISTS CONTAINER ##############à
  containername="clistNCsign";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCsign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistNCsign);
  last++;


  containername="clistNCback";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCback = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistNCback);
  last++;

  containername="clistNCfromB";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCfromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistNCfromB);
  last++;


  containername="clistNCfromDstar";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCfromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistNCfromDstar);
  last++;


  containername="clistNCother";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCother = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistNCother);
  last++;


  //######### LOOSE CUTS TLISTS CONTAINER #############
  containername="clistLSCsign";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCsign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCsign);
  last++;


  containername="clistLSCback";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCback = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCback);
  last++;

  containername="clistLSCfromB";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCfromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCfromB);
  last++;


  containername="clistLSCfromDstar";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCfromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCfromDstar);
  last++;


  containername="clistLSCother";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCother = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCother);
  last++;



  //######### TIGHT CUTS TLISTS CONTAINER #############
    containername="clistTGHCsign";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCsign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCsign);
  last++;


  containername="clistTGHCback";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCback = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCback);
  last++;

  containername="clistTGHCfromB";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCfromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCfromB);
  last++;


  containername="clistTGHCfromDstar";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCfromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCfromDstar);
  last++;


  containername="clistTGHCother";
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCother = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outfile.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCother);
  


  return hfTask;
}
