AliAnalysisTaskSEXicZero2XiPifromAODtracks *AddTaskXicZero2XiPifromAODtracks(TString finname="",
									     TString outputFileName="",
										 Bool_t theMCon=kFALSE,
                     Bool_t anaOmegacZero=kTRUE,
										 Bool_t writeVariableTree=kTRUE,
										 Bool_t reconstructPrimVert=kFALSE,
									         Bool_t fillOnlySig = kFALSE,
										 Bool_t fillOnlyBkg = kFALSE
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
      cout << "Input file not found : check your cut object" << endl;
      return 0x0;
      //      AliFatal("Input file not found : check your cut object");
    }
  }
  
  AliRDHFCutsXicZerotoXiPifromAODtracks* RDHFCutsXic2XiPianal= new AliRDHFCutsXicZerotoXiPifromAODtracks();
  if (stdcuts) RDHFCutsXic2XiPianal->SetStandardCutsPP2010();
  else RDHFCutsXic2XiPianal = (AliRDHFCutsXicZerotoXiPifromAODtracks*)filecuts->Get("XicZeroAnalysisCuts");
  RDHFCutsXic2XiPianal->SetName("XicZeroAnalysisCuts");
  RDHFCutsXic2XiPianal->SetMinPtCandidate(2.);
  RDHFCutsXic2XiPianal->SetMaxPtCandidate(10000.);

    
  // mm let's see if everything is ok
  if (!RDHFCutsXic2XiPianal) {
    cout << "Specific AliRDHFCutsXic2XiPianal not found\n";
    return 0x0;
  }


  //CREATE THE TASK
  
  printf("CREATE TASK\n");
  AliAnalysisTaskSEXicZero2XiPifromAODtracks *task = new AliAnalysisTaskSEXicZero2XiPifromAODtracks("AliAnalysisTaskSEXicZero2XiPifromAODtracks", RDHFCutsXic2XiPianal, writeVariableTree, anaOmegacZero);
  task->SetMC(theMCon);
  task->SetFillSignalOnly(fillOnlySig);
  task->SetFillBkgOnly(fillOnlyBkg);
  if ((!theMCon && fillOnlySig) || (!theMCon && fillOnlyBkg) || (fillOnlySig && fillOnlyBkg)) {cout<< "Wrong settings for MC" << endl; return 0x0;}
  task->SetDebugLevel(1);
  task->SetReconstructPrimVert(reconstructPrimVert);
  
  mgr->AddTask(task);
  
  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if (!anaOmegacZero) outputfile += ":PWG3_D2H_XicZero2XiPi_";
  else outputfile += ":PWG3_D2H_OmegacZero2OmegaPi_";
  outputfile += outputFileName.Data();
  
  // Connect output of an existing task to a data container
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  
  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("chist%s", outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);

  AliAnalysisDataContainer *coutputXic2 = NULL;
  if (!anaOmegacZero)
    coutputXic2 = mgr->CreateContainer(Form("XicZero2XiPiCuts%s", outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  else
    coutputXic2 = mgr->CreateContainer(Form("OmegacZero2OmegaPiCuts%s", outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  mgr->ConnectOutput(task,2,coutputXic2);

  if (writeVariableTree) {
    AliAnalysisDataContainer *coutputXic3 = NULL;
    if (!anaOmegacZero)
      coutputXic3 = mgr->CreateContainer(Form("XicZerovariables%s", outputFileName.Data()),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
    else
      coutputXic3 = mgr->CreateContainer(Form("OmegacZerovariables%s", outputFileName.Data()),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    mgr->ConnectOutput(task,3,coutputXic3);
  }else{
    AliAnalysisDataContainer *coutputXic3 = NULL;
    if (!anaOmegacZero)
      coutputXic3 = mgr->CreateContainer(Form("XicZeroAll%s", outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
    else
      coutputXic3 = mgr->CreateContainer(Form("OmegacZeroAll%s", outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
    mgr->ConnectOutput(task,3,coutputXic3);
  }
  
  return task;

}
