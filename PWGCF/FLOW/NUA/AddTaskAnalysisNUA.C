//_________________________________________________________//
AliAnalysisTaskNUA *AddTaskAnalysisNUA(Double_t vertexZ=10.,
						UInt_t triggerSelectionString = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
						Int_t gFilterBit = 768,
						Double_t etaMin=-0.8, Double_t etaMax=0.8,
						TString fileNameBase="AnalysisResults") {
  // Creates an analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  TString centralityName("");

  static const Int_t nCentralities = 12;
  Double_t gCentrality[nCentralities];

  gCentrality[0]=0;
  gCentrality[1]=5;
  for(Int_t i = 2; i <= 11; i++)
  gCentrality[i] = 0 + (i-1)*10;
     
  //===========================================================================

  AliAnalysisTaskNUA *task[nCentralities];
  AliAnalysisDataContainer *coutQA[nCentralities];
  AliAnalysisDataContainer *coutResults[nCentralities];
  TString suffixName[nCentralities];

  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {

    suffixName[iCentralityBin] = "Centrality";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin];
    suffixName[iCentralityBin] += "To";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAnalysisAOD", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskAnalysisAOD", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Create the task, add it to manager and configure it.
  //===========================================================================
  //AliAnalysisTaskAODPtSpectra *task = new AliAnalysisTaskAODPtSpectra("TaskAODPtSpectra");
  // offline trigger selection (AliVEvent.h)
   
  task[iCentralityBin] = new AliAnalysisTaskNUA(Form("Task_%s",suffixName[iCentralityBin].Data()));
  task[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);

 //Centrality estimator

  task[iCentralityBin]->SetCentralityEstimator("V0M");
  
  // vertex cut (x,y,z)
  task[iCentralityBin]->SetVertexDiamond(30.,30.,vertexZ);
  
  // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
  task[iCentralityBin]->SetAODtrackCutBit(gFilterBit);
  task[iCentralityBin]->SetKinematicsCutsAOD(etaMin,etaMax);
  task[iCentralityBin]->SetCentralityPercentileRange(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1]);
  task[iCentralityBin]->SetMinPt(0);
  task[iCentralityBin]->SetMaxPt(10.0);

  mgr->AddTask(task[iCentralityBin]);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //===========================================================================
  mgr->ConnectInput(task[iCentralityBin], 0, mgr->GetCommonInputContainer());

  coutQA[iCentralityBin] = mgr->CreateContainer(Form("listQA_%s", suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectOutput(task[iCentralityBin], 1, coutQA[iCentralityBin]);

  coutResults[iCentralityBin] = mgr->CreateContainer(Form("listResults_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectOutput(task[iCentralityBin], 2, coutResults[iCentralityBin]);
   
  //return task;
}

}
