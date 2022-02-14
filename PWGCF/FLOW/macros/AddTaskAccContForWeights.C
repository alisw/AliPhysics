//_________________________________________________________//    
    
AliAnalysisTaskAccCont *AddTaskAccContForWeights(Double_t vertexZ=10.,
				       UInt_t triggerSelectionString = AliVEvent::kINT7, //AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
				       Int_t gFilterBit = 768, Int_t nsigma = 3,
				       Double_t ptMin=0, Double_t ptMax=10,
				       Double_t etaMin=-0.8, Double_t etaMax=0.8,
				       Bool_t DCAext = kFALSE,
				       Bool_t PID = kFALSE,
				       Bool_t UseRapidity = kFALSE,
				       TString centEstim = "V0M",
				       AliAnalysisTaskAccCont::kSystem systemType = AliAnalysisTaskAccCont::kPbPb,
				       AliAnalysisTaskAccCont::kCentralityBinning nCenBins = AliAnalysisTaskAccCont::kBins, //AliAnalysisTaskAccCont::kFull,
				       AliPID::EParticleType particleType = AliPID::kMuon,
				       TString fileNameBase="AnalysisResults",
		                       TString MCorData = "data") {
  // Creates an analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  TString centralityName("");

  static const Int_t nCentralitiesMax = 9;
  Double_t gCentrality[nCentralitiesMax];

  if(nCenBins == AliAnalysisTaskAccCont::kFull){

  gCentrality[0]=0;
  gCentrality[1]=80;

  }

  else if (nCenBins == AliAnalysisTaskAccCont::kBins){
  
  for(Int_t i = 0; i < 9; i++){
  gCentrality[i] = 0 + i*10;
  }
  
  }

  else if(nCenBins == AliAnalysisTaskAccCont::kMCgen){

  gCentrality[0]=0;
  gCentrality[1]=100;

  }
  //===========================================================================


  Int_t nCentralities = 2;
  if(nCenBins == AliAnalysisTaskAccCont::kFull || nCenBins == AliAnalysisTaskAccCont::kMCgen) nCentralities = 2;
  else if(nCenBins == AliAnalysisTaskAccCont::kBins) nCentralities = nCentralitiesMax;

  AliAnalysisTaskAccCont *task[nCentralitiesMax];
  AliAnalysisDataContainer *coutQA[nCentralitiesMax];
  AliAnalysisDataContainer *coutResults[nCentralitiesMax];
  TString suffixName[nCentralitiesMax];

  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
	std::cout<<"================================================> Start centrality Bin loop: "<<iCentralityBin<<std::endl;
	
    suffixName[iCentralityBin] = "Centrality";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin];
    suffixName[iCentralityBin] += "To";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];
    suffixName[iCentralityBin] += "_fb";
    suffixName[iCentralityBin] += gFilterBit;
    suffixName[iCentralityBin] += "_CentEst_";
    suffixName[iCentralityBin] += centEstim.Data();
	
    if(PID){
    suffixName[iCentralityBin] += "_PID_";
    suffixName[iCentralityBin] += AliPID::ParticleShortName(particleType);
    }
    
    if(UseRapidity)
    suffixName[iCentralityBin] += "_RapidityUsed"; 

    if(nCenBins == AliAnalysisTaskAccCont::kFull)
    suffixName[iCentralityBin] += "_full";
    else if(nCenBins == AliAnalysisTaskAccCont::kBins)
    suffixName[iCentralityBin] += "_bins";
    else if(nCenBins == AliAnalysisTaskAccCont::kMCgen)
    suffixName[iCentralityBin] += "_MCgen";
	std::cout<<"suffixName["<<iCentralityBin<<"] ==== "<<suffixName[iCentralityBin]<<std::endl;

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
   
  task[iCentralityBin] = new AliAnalysisTaskAccCont(Form("Task_%s",suffixName[iCentralityBin].Data()));
  task[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
  
  //Centrality estimator
  task[iCentralityBin]->SetCentralityEstimator(centEstim.Data());
  
  //vertex cut (x,y,z)
  task[iCentralityBin]->SetVertexDiamond(3.,3.,vertexZ);
  
  //pt and eta cut (pt_min, pt_max, eta_min, eta_max)
  task[iCentralityBin]->SetAODtrackCutBit(gFilterBit);
  task[iCentralityBin]->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
  task[iCentralityBin]->SetCentralityPercentileRange(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1]);
  
  if(systemType == AliAnalysisTaskAccCont::kPbPb){
  if (MCorData == "data"){
  task[iCentralityBin]->CheckPileUp();
  task[iCentralityBin]->UsePileUpCutsPbPb();
  }
  else if (MCorData == "MCrec")
  task[iCentralityBin]->SetMCRec();
  }
  else if(systemType == AliAnalysisTaskAccCont::kpPb){
  if (MCorData == "data"){
  task[iCentralityBin]->CheckPileUp();
  task[iCentralityBin]->UsePileUpCutspPb();
  }
  else if (MCorData == "MCrec")
  task[iCentralityBin]->SetMCRec();
  }

  if(UseRapidity)
  task[iCentralityBin]->SetUseRapidity(); 
 
  if(DCAext)
  task[iCentralityBin]->USEextendedDCA();
  
  
  //PID
  if(PID){
  if(systemType == AliAnalysisTaskAccCont::kpPb)
  task[iCentralityBin]->SetUseNSigmaPIDNewTrial();
  else
  task[iCentralityBin]->UsePID();
  //task[iCentralityBin]->SetNSigmaPID(nsigma);
  task[iCentralityBin]->setParticleType(particleType);
  mgr->AddTask(task[iCentralityBin]);  
  }
  
  else if(!PID)
  mgr->AddTask(task[iCentralityBin]);
  

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //===========================================================================
  mgr->ConnectInput(task[iCentralityBin], 0, mgr->GetCommonInputContainer());

  coutQA[iCentralityBin] = mgr->CreateContainer(Form("listQA_%s", suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectOutput(task[iCentralityBin], 1, coutQA[iCentralityBin]);

  coutResults[iCentralityBin] = mgr->CreateContainer(Form("listResults_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectOutput(task[iCentralityBin], 2, coutResults[iCentralityBin]);

  //return task[iCentralityBin];
  
  }
  return task[0];
}
