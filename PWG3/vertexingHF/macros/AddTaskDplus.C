AliAnalysisTaskSEDplus *AddTaskDplus(Bool_t storeNtuple=kFALSE,
				     Bool_t readMC=kFALSE,
				     TString filename="DplustoKpipiCuts.root")
{
  //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for D+ candidates 

  //Invariant mass histogram and                                                 
  // association with MC truth (using MC info in AOD)                                                                                   
  //  R. Bala, bala@to.infn.it                                                                                                                                  
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDplus", "No analysis manager to connect to.");
  }

  Bool_t stdcuts=kFALSE;
  TFile* filecuts=new TFile(filename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: using standard cuts"<<endl;
    stdcuts=kTRUE;
  }
  
  
  //Analysis Task

  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  if(stdcuts) analysiscuts->SetStandardCutsPP2010();
  else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get("AnalysisCuts");

  AliRDHFCutsDplustoKpipi* prodcuts=new AliRDHFCutsDplustoKpipi();
  if(stdcuts) prodcuts->SetStandardCutsPP2010();
  else prodcuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get("ProdCuts");
  
  //AliRDHFCutsDplustoKpipi *prodcuts = (AliRDHFCutsDplustoKpipi*)fileCuts->Get("ProdCuts");
  //AliRDHFCutsDplustoKpipi *analysiscuts = (AliRDHFCutsDplustoKpipi*)fileCuts->Get("AnalysisCuts");

  
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",analysiscuts,prodcuts,storeNtuple);
  dplusTask->SetReadMC(readMC);
  dplusTask->SetDoLikeSign(kFALSE);
  //  dplusTask->SetUseTPCpid(kTRUE);
  //dplusTask->SetUseTOFpid(kTRUE);
  dplusTask->SetDebugLevel(0);
  dplusTask->SetMassLimits(0.2);
  dplusTask->SetUseBit(kTRUE);
  mgr->AddTask(dplusTask);
  
  // Create containers for input/output 
  
  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer("cinputDplus",TChain::Class(),
							       AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDplus";
  
  AliAnalysisDataContainer *coutputDplusCuts = mgr->CreateContainer("coutputDplusCuts",TList::Class(),
								    AliAnalysisManager::kOutputContainer,
								    outputfile.Data());
  
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer("coutputDplus",TList::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputDplusNorm = mgr->CreateContainer("coutputDplusNorm",AliNormalizationCounter::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  
  if(storeNtuple){
    AliAnalysisDataContainer *coutputDplus2 = mgr->CreateContainer("coutputDplus2",TNtuple::Class(),
								   AliAnalysisManager::kOutputContainer,
								   "InvMassDplus_nt1.root");
    
    coutputDplus2->SetSpecialOutput();
  }
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dplusTask,1,coutputDplus);
  
  mgr->ConnectOutput(dplusTask,2,coutputDplusCuts);

  mgr->ConnectOutput(dplusTask,3,coutputDplusNorm);  
  if(storeNtuple){
    mgr->ConnectOutput(dplusTask,4,coutputDplus2);
  }
  return dplusTask;
}
