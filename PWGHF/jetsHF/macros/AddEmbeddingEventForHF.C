AliEmbeddingEventForHFTask *AddEmbeddingEventForHF(
  TString BaseTrackArrName = "",
  TString aodTreeName = "aodTree",
  TString aodTrackArrName = "tracks",
  TString gridDir = "",
  TString passAndDataType = "",
  TString dataFile = "",
  Int_t maxFile = 10,
  Bool_t isSimulation = kFALSE, //if the sample to be embedded is simulation, not the base sample!
  Double_t minCent = 0,
  Double_t maxCent = 10,
  UInt_t mask = AliVEvent::kCentral,
  TString taskname = "",
  TString cutFile = "DStartoKpipiCuts.root",
  AliEmbeddingEventForHFTask::ACandidateType cand = AliEmbeddingEventForHFTask::kDstartoKpipi 
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddEmbeddingEventForHF", "No analysis manager to connect to.");
    return NULL;
  } 

    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskJetEmbeddingFromAOD", "This task requires an input event handler");
        return NULL;
    }

    Bool_t useStdC = kFALSE;
    TFile* filecuts = TFile::Open(cutFile);
    if (!filecuts || (filecuts && !filecuts->IsOpen())) {
    ::Warning("AddTaskSEDmesonsFilterCJ1", "Input file not found: use std cuts");
    useStdC = kTRUE;
  }


  TString candname("DStar"); 
  if (cand==0) candname="D0";

  AliRDHFCuts *analysiscuts = 0x0;
  switch (cand) {
  case 0 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
    break;
  case 1 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    analysiscuts->SetName("DStartoKpipiCuts");
    break;
  }
  
  if (!analysiscuts) {
    std::cout<<"Specific AliRDHFCuts not found"<<std::endl;
    return;
  }
    
    AliEmbeddingEventForHFTask* taskEmb = mgr->GetTask(taskname.Data());
    
    if (taskEmb) {
    ::Info("AddEmbeddingEventForHF", Form("Task %s already exist, continue",taskname.Data()));
    return taskEmb;
  }
  else
  {
    ::Info("AddEmbeddingEventForHF", "Creating the task");
    
    taskEmb = new AliEmbeddingEventForHFTask(taskname.Data(),analysiscuts,cand);
    
    taskEmb->SetGridDir(gridDir);
    taskEmb->SetPassAndDataType(passAndDataType);
    taskEmb->SetDataFile(dataFile);
    taskEmb->SetMaxFileNumber(maxFile);
    taskEmb->SetGridDir(isSimulation);
    
    taskEmb->SetAODTreeName(aodTreeName);
    taskEmb->SetAODTracksName(aodTrackArrName);
    taskEmb->SetCentralityRange(minCent, maxCent);
    taskEmb->SetTriggerMask(mask);
    
    mgr->AddTask(taskEmb);
    
  }
    
    // Create and connect containers for input/output
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":Embedding";
    
    TString nameContainer0 = "histograms";
    nameContainer0 += taskname;
   
    // ------ input data ------
    AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
    
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(nameContainer0, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());

    
    
    mgr->ConnectInput(taskEmb, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskEmb,1,coutput1);
    
    ::Info("AddEmbeddingEventForHF", "Input and Output connected to the manager");
    return taskEmb;
    
}
