AliAnalysisTaskSEDs *AddTaskDs(Int_t storeNtuple=0,Bool_t readMC=kFALSE,
				     TString filename="DstoKKpiCuts.root")
{
  //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for Ds candidates 

  //Invariant mass histogram and                                                 
  // association with MC truth (using MC info in AOD)                                                                                   
  // Origin: R. Bala, bala@to.infn.it         
  // Modified for Ds meson: G.M. Innocenti innocent@to.infn.it
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDs", "No analysis manager to connect to.");
  }


  TFile* filecuts=new TFile(filename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Error: Input file not found!"<<endl;
    return 0;
  }
  
  
  //Analysis Task

  
  AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
  analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get("AnalysisCuts");

  AliRDHFCutsDstoKKpi* prodcuts=new AliRDHFCutsDstoKKpi();
  prodcuts = (AliRDHFCutsDstoKKpi*)filecuts->Get("ProdCuts");
  
  //AliRDHFCutsDstoKKpi *prodcuts = (AliRDHFCutsDstoKKpi*)fileCuts->Get("ProdCuts");
  //AliRDHFCutsDstoKKpi *analysiscuts = (AliRDHFCutsDstoKKpi*)fileCuts->Get("AnalysisCuts");

  AliAnalysisTaskSEDs *dsTask = new AliAnalysisTaskSEDs("DsAnalysis",prodcuts,analysiscuts,storeNtuple);
  dsTask->SetReadMC(readMC);
  //dsTask->SetDoLikeSign(kTRUE);
  //  dsTask->SetUseTPCpid(kTRUE);
  //dsTask->SetUseTOFpid(kTRUE);
  dsTask->SetDebugLevel(10);
  //dsTask->SetMassLimits(0.2);
  mgr->AddTask(dsTask);
  
  // Create containers for input/output 
  
  AliAnalysisDataContainer *cinputDs = mgr->CreateContainer("cinputDs",TChain::Class(),
							       AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDs";
  
  AliAnalysisDataContainer *coutputDsCuts = mgr->CreateContainer("coutputDsCuts",TList::Class(),
								    AliAnalysisManager::kOutputContainer,
								    outputfile.Data());
  
  AliAnalysisDataContainer *coutputDs = mgr->CreateContainer("coutputDs",TList::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputDsNorm = mgr->CreateContainer("coutputDsNorm",AliNormalizationCounter::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  
   if(storeNtuple){
    AliAnalysisDataContainer *coutputDs2 = mgr->CreateContainer("coutputDs2",TNtuple::Class(),
								   AliAnalysisManager::kOutputContainer,
								   "InvMassDs_nt1.root");
    
    coutputDs2->SetSpecialOutput();
  }
  
  mgr->ConnectInput(dsTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dsTask,1,coutputDs);
  
  mgr->ConnectOutput(dsTask,2,coutputDsCuts);

  mgr->ConnectOutput(dsTask,3,coutputDsNorm);  
  
  if(storeNtuple){
    mgr->ConnectOutput(dsTask,4,coutputDs2);
  }
  
  return dsTask;
}
