AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks *AddTaskXicPlus2XiPiPifromAODtracks(TString finname="",
										 TString outputFileName="",
										 Bool_t theMCon=kFALSE,
										 Bool_t writeVariableTree=kTRUE,
										 Bool_t reconstructPrimVert=kFALSE,
										 Bool_t fillOnlySig = kFALSE,
										 Bool_t fillOnlyBkg = kFALSE,
										 Bool_t fillSparse = kFALSE,
										 Bool_t HMTrigOn = kTRUE, //jcho, for HM analysis
										 Bool_t EvtInfo = kFALSE  //jcho, to write the event tree
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
  
  AliRDHFCutsXicPlustoXiPiPifromAODtracks* RDHFCutsXic2PiPianal = new AliRDHFCutsXicPlustoXiPiPifromAODtracks();
  if (stdcuts) RDHFCutsXic2PiPianal->SetStandardCutsPP2010();
  else RDHFCutsXic2PiPianal = (AliRDHFCutsXicPlustoXiPiPifromAODtracks*)filecuts->Get("XicPlusAnalysisCuts");
  RDHFCutsXic2PiPianal->SetName("XicPlusAnalysisCuts");
  RDHFCutsXic2PiPianal->SetMinPtCandidate(2.); 
  RDHFCutsXic2PiPianal->SetMaxPtCandidate(10000.);

    
  // mm let's see if everything is ok
  if (!RDHFCutsXic2PiPianal) {
    cout << "Specific AliRDHFCutsXic2PiPianal not found\n";
    return 0x0;
  }


  //CREATE THE TASK
  if(writeVariableTree && fillSparse) {
    cout << "You are enabling the filling of both Tree and THnSparse\n";
    cout << "Please choose only one of the two options.\n";
    return 0x0;
  }
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks *task = new AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks("AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks",RDHFCutsXic2PiPianal,writeVariableTree,fillSparse,HMTrigOn,EvtInfo);  
  task->SetMC(theMCon);
  task->SetFillSignalOnly(fillOnlySig);
  task->SetFillBkgOnly(fillOnlyBkg);
  if ((!theMCon && fillOnlySig) || (!theMCon && fillOnlyBkg) || (fillOnlySig && fillOnlyBkg)) {cout<< "Wrong settings for MC" << endl; return 0x0;}
  task->SetDebugLevel(1);
  task->SetReconstructPrimVert(reconstructPrimVert);
  
  //mgr->AddTask(task);  
  task->UseTrig_kINT7();	//jcho
  if (HMTrigOn){
	task->UseTrig_kHMV0();	//jcho, HM Trigger
	//task->UseTrig_kHMSPD();  //jcho, HM Trigger
	}

  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_XicPlus2XiPiPi_";
  outputfile += outputFileName.Data();
  if(HMTrigOn) outputfile += "_HM";

  mgr->AddTask(task); //jcho
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  
  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("chist%s",outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutputXic2 = mgr->CreateContainer(Form("XicPlus2XiPiPiCuts%s",outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  mgr->ConnectOutput(task,2,coutputXic2);
  if (writeVariableTree) {
    AliAnalysisDataContainer *coutputXic3 = mgr->CreateContainer(Form("XicPlusvariables%s",outputFileName.Data()),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
    mgr->ConnectOutput(task,3,coutputXic3);
  }else{
    AliAnalysisDataContainer *coutputXic3 = mgr->CreateContainer(Form("XicPlusAll%s",outputFileName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
    mgr->ConnectOutput(task,3,coutputXic3);
  }
   AliAnalysisDataContainer *coutputXic4 = mgr->CreateContainer(Form("XicPlus2XiPiPiNorm%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // normalization counter
  mgr->ConnectOutput(task,4,coutputXic4);

  // For HM analysis -------jcho
  AliAnalysisDataContainer *coutputXic5 = mgr->CreateContainer(Form("MB_0to100%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
  mgr->ConnectOutput(task,5,coutputXic5);
  AliAnalysisDataContainer *coutputXic6 = mgr->CreateContainer(Form("MB_0point1to30%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
  mgr->ConnectOutput(task,6,coutputXic6);
  AliAnalysisDataContainer *coutputXic7 = mgr->CreateContainer(Form("MB_30to100%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,7,coutputXic7);
  AliAnalysisDataContainer *coutputXic8 = mgr->CreateContainer(Form("HMV0_0to0point1%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,8,coutputXic8); 
  AliAnalysisDataContainer *coutputXic9 = mgr->CreateContainer(Form("HMV0_0to100%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,9,coutputXic9);

  //
  AliAnalysisDataContainer *coutputXic10 = mgr->CreateContainer(Form("MB_0to100_INEL%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,10,coutputXic10);
  AliAnalysisDataContainer *coutputXic11 = mgr->CreateContainer(Form("MB_0point1to30_INEL%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,11,coutputXic11);
  AliAnalysisDataContainer *coutputXic12 = mgr->CreateContainer(Form("MB_30to100_INEL%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,12,coutputXic12);
  AliAnalysisDataContainer *coutputXic13 = mgr->CreateContainer(Form("HMV0_0to100_INEL%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,13,coutputXic13);
  AliAnalysisDataContainer *coutputXic14 = mgr->CreateContainer(Form("HMV0_0to0point1_INEL%s",outputFileName.Data()),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,14,coutputXic14);

  if(theMCon){ 
  AliAnalysisDataContainer *coutputXic15 = mgr->CreateContainer(Form("MCGenVariables%s",outputFileName.Data()),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(task,15,coutputXic15);
  }

  if(EvtInfo) {
  AliAnalysisDataContainer *coutputXic16 = mgr->CreateContainer(Form("XicPlusEventvariables%s",outputFileName.Data()),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());    //jcho, Event variables tree
  mgr->ConnectOutput(task,16,coutputXic16);
  }

  return task;

}
