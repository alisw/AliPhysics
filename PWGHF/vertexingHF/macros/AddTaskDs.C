AliAnalysisTaskSEDs *AddTaskDs(Int_t system=0/*0=pp,1=PbPb*/,
                     Int_t storeNtuple=0,Bool_t readMC=kFALSE,
				     TString filename="")
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
  
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
    filecuts=TFile::Open(filename.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	  AliFatal("Cut object not found: analysis will not start!\n");
    }
    else printf("Cut object correctly found\n");
  }
  
  //Analysis Task

  AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
  
  if(stdcuts) {
    if(system==0) {
      printf("Cut object not found: standard pp cut object used\n");
      analysiscuts->SetStandardCutsPP2010();
    }
    else AliFatal("Standard cut object not available for PbPb: analysis will not start!\n");
  }
  else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get("AnalysisCuts");
  
  AliAnalysisTaskSEDs *dsTask = new AliAnalysisTaskSEDs("DsAnalysis",analysiscuts,storeNtuple);

  dsTask->SetReadMC(readMC);
  //dsTask->SetDoLikeSign(kTRUE);
  //  dsTask->SetUseTPCpid(kTRUE);
  //dsTask->SetUseTOFpid(kTRUE);
  dsTask->SetDebugLevel(10);
  dsTask->SetUseSelectionBit(kTRUE);
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
								   outputfile.Data());
    
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
