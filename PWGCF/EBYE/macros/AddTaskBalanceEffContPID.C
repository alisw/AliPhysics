AliAnalysisTaskEffContPIDBF *AddTaskBalanceEffContPID(TString  centralityEstimator="V0M",
						 Double_t centrMin=0.,
						 Double_t centrMax=80.,
						 Double_t vertexZ=10.,
                                                 Double_t ptmin = 0.2,
                                                 Double_t ptmax = 2.0,
                                                 Double_t pttpcmax = 0.6,
                                                 Double_t nsigmapid = 3.0,
                                                 Int_t pidType = 1,
						 Int_t AODfilterBit = 128,
						 Bool_t bUseElectronRejection = kFALSE,
						 TString fileNameBase="AnalysisResults",
                                                 TString dirNameExtra = ""
						 ) {

  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");

  centralityName+=Form("%.0f-%.0f_%.0f",centrMin,centrMax);
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTriggeredBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTriggeredBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";

 
  // Create the task, add it to manager and configure it.
  //===========================================================================


  AliAnalysisTaskEffContPIDBF *taskEffContPIDBF = new AliAnalysisTaskEffContPIDBF("TaskEffContPIDBF");

  // centrality
  if(centralityEstimator) {
    taskEffContPIDBF->UseCentrality();
    taskEffContPIDBF->SetCentralityEstimator(centralityEstimator);
    taskEffContPIDBF->SetCentralityPercentileRange(centrMin,centrMax);
  }
  taskEffContPIDBF->SelectCollisionCandidates(AliVEvent::kMB);

  // vertex
  taskEffContPIDBF->SetVertexDiamond(.3,.3,vertexZ);

  //analysis kinematic cuts
  taskEffContPIDBF->SetMinPt(ptmin);
  taskEffContPIDBF->SetMaxPt(ptmax); //5.0
   taskEffContPIDBF->SetTPCPtMax(pttpcmax); // For TPC Max
  //taskEffContPIDBF->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
  //taskEffContPIDBF->SetPtRange(0.1, 20.0, 100);  //acceptance cuts //5.0,49
  taskEffContPIDBF->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
  taskEffContPIDBF->SetPtRange(0.2, 20.0, 100);  //acceptance cuts //5.0,49
  taskEffContPIDBF->SetNSigmaCut(nsigmapid); // NSigma Cut;
  taskEffContPIDBF->SetPIDType(pidType); // KNSigmaTPCTOF == 1, KNSigmaTPC==0, KNSigmaTOF == 2
  


  // electron rejection
    if(bUseElectronRejection){
      taskEffContPIDBF->SetElectronOnlyRejection(3.); // no other particle in nsigma (this is what we use standard in BF code)
    }

  //AODs  
  taskEffContPIDBF->SetAODtrackCutBit(AODfilterBit);
  mgr->AddTask(taskEffContPIDBF);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionEffContPIDAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_%s_Bit%d_%s%s",centralityName.Data(),AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutEffContBF = mgr->CreateContainer(Form("listEffContBF_%s_Bit%d_%s%s",centralityName.Data(),AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectInput(taskEffContPIDBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEffContPIDBF, 1, coutQA);
  mgr->ConnectOutput(taskEffContPIDBF, 2, coutEffContBF);
  
  return taskEffContPIDBF;
}
