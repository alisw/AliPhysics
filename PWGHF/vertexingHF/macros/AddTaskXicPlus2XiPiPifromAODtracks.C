AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks *AddTaskXicPlus2XiPiPifromAODtracks(TString finname="",
										 Bool_t theMCon=kFALSE,
										 Bool_t writeVariableTree=kFALSE,
										 Int_t nTour=0
										 )

{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLc2V0YW", "No analysis manager to connect to.");
    return NULL;
  }  
  
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( finname.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
    filecuts=TFile::Open(finname.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      AliFatal("Input file not found : check your cut object");
    }
  }
  AliRDHFCutsXicPlustoXiPiPifromAODtracks* RDHFCutsXic2PiPiprod = new AliRDHFCutsXicPlustoXiPiPifromAODtracks();
  if (stdcuts) RDHFCutsXic2PiPiprod->SetStandardCutsPP2010();
  else RDHFCutsXic2PiPiprod = (AliRDHFCutsXicPlustoXiPiPifromAODtracks*)filecuts->Get("XicPlusProductionCuts");
  RDHFCutsXic2PiPiprod->SetName("XicPlusProductionCuts");
  RDHFCutsXic2PiPiprod->SetMinPtCandidate(2.);
  RDHFCutsXic2PiPiprod->SetMaxPtCandidate(10000.);

  if (!RDHFCutsXic2PiPiprod) {
    cout << "Specific AliRDHFCutsXic2PiPiprod not found\n";
    return;
  }
  
  AliRDHFCutsXicPlustoXiPiPifromAODtracks* RDHFCutsXic2PiPianal = new AliRDHFCutsXicPlustoXiPiPifromAODtracks();
  if (stdcuts) RDHFCutsXic2PiPianal->SetStandardCutsPP2010();
  else RDHFCutsXic2PiPianal = (AliRDHFCutsXicPlustoXiPiPifromAODtracks*)filecuts->Get("XicPlusAnalysisCuts");
  RDHFCutsXic2PiPianal->SetName("XicPlusAnalysisCuts");
  RDHFCutsXic2PiPianal->SetMinPtCandidate(2.);
  RDHFCutsXic2PiPianal->SetMaxPtCandidate(10000.);

    
  // mm let's see if everything is ok
  if (!RDHFCutsXic2PiPianal) {
    cout << "Specific AliRDHFCutsXic2PiPianal not found\n";
    return;
  }


  //CREATE THE TASK
  
  printf("CREATE TASK\n");
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks *task = new AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks("AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks",RDHFCutsXic2PiPiprod,RDHFCutsXic2PiPianal,writeVariableTree);
  task->SetIspp(kFALSE);
  task->SetIspA(kFALSE);
  task->SetIsAA(kTRUE);
  task->SetMC(theMCon);
  task->SetDebugLevel(1);
  
  mgr->AddTask(task);
  
  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_XicPlus2XiPiPi_";
  outputfile += nTour;
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  
  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("chist%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutputXic2 = mgr->CreateContainer(Form("XicPlus2XiPiPiCuts%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  mgr->ConnectOutput(task,2,coutputXic2);
  if (writeVariableTree) {
    AliAnalysisDataContainer *coutputXic3 = mgr->CreateContainer(Form("XicPlusvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
    mgr->ConnectOutput(task,3,coutputXic3);
  }else{
    AliAnalysisDataContainer *coutputXic3 = mgr->CreateContainer(Form("XicPlusAll%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
    mgr->ConnectOutput(task,3,coutputXic3);
  }
  
  return task;

}
