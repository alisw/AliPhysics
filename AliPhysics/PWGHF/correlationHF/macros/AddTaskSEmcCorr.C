AliAnalysisTaskSEmcCorr* AddTaskSEmcCorr(TString fileout="d0D0.root",TString containerprefix="c",Bool_t readmc=kTRUE,Bool_t doHH=kFALSE,TString genTitle="")
{  
  //
  // andrea.rossi@cern.ch
  //
  //==========================================================================

  Int_t last=0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCheckNprongs", "No analysis manager to connect to.");
    return NULL;
  }   
  
  TString str,containername;
  if(fileout=="standard"){
    fileout=AliAnalysisManager::GetCommonFileName();
    fileout+=":PWG3_HFCK_";
    fileout+="TestNprongs";
    if(containerprefix!="c")fileout+=containerprefix;
    str="TestNprongs";
  }
  else {
    str=fileout;
    str.ReplaceAll(".root","");
  }
  str.Prepend("_");

  AliAnalysisTaskSEmcCorr *hfTask = new AliAnalysisTaskSEmcCorr("AliAnalysisTaskSEmcCorr");
  hfTask->SetReadMC(readmc);
  hfTask->DoHadronHadron(doHH);
  hfTask->SetGeneratorToBeChecked(genTitle);
  mgr->AddTask(hfTask);
 
 
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput =   mgr->GetCommonInputContainer();
  //mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,cinput);
  

  //Now container for general properties histograms
  containername="outputNentries";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,1,coutputNentries);

  containername="outputNprongsD0";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNprongsD0 = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,2,coutputNprongsD0);


  containername="outputNprongsD0chargedOnly";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNprongsD0chargedOnly = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,3,coutputNprongsD0chargedOnly);

  containername="outputNprongsD0chargedRef";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNprongsD0chargedRef = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,4,coutputNprongsD0chargedRef);


  containername="outputMCcorrel";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputMCcorrel = mgr->CreateContainer(containername.Data(),THnSparseF::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,5,coutputMCcorrel);


  containername="outputMCcorrelTrig";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputMCcorrelTrig = mgr->CreateContainer(containername.Data(),THnSparseF::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,6,coutputMCcorrelTrig);



  containername="outputMChadroncorrel";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputMChadroncorrel = mgr->CreateContainer(containername.Data(),THnSparseF::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,7,coutputMChadroncorrel);


  containername="outputMChadroncorrelTrig";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputMChadroncorrelTrig = mgr->CreateContainer(containername.Data(),THnSparseF::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,8,coutputMChadroncorrelTrig);


  containername="outputDzeroDecay";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputDzeroDecay = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,9,coutputDzeroDecay);


  containername="outputDplusDecay";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputDplusDecay = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,10,coutputDplusDecay);


  containername="outputLambdaCDecay";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputLambdaCDecay = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,11,coutputLambdaCDecay);



  containername="outputAllBDecay";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputAllBDecay = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,12,coutputAllBDecay);
  return hfTask;
}
