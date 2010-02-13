//DEFINITION OF A FEW CONSTANTS
const Int_t    mintrackrefsTPC = 0 ;
const Int_t    mintrackrefsITS = 3 ;
const Int_t    PDG = 421; 
const Int_t    minclustersTPC = 0 ;
const Int_t    minITSClusters = 4;
// ANALYSIS TYPE D*+ or D*-
const Bool_t computeD0 = kFALSE;
const Bool_t topologicalCut = kFALSE;
//----------------------------------------------------

AliAnalysisTaskSEDStarJets *AddTaskDStarJets(Bool_t theMCon=kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDStarJets", "No analysis manager to connect to.");
    return NULL;
  }  
  
  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliAnalysisTaskSEDStarJets *task = new AliAnalysisTaskSEDStarJets("AliAnalysisTaskSEDStarJets");
  task->SetMinITSClusters(minITSClusters);
  task->SetAnalType(computeD0);
  task->SetMC(theMCon);
  task->SetCutType(topologicalCut);
  // Create and connect containers for input/output
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DStarJet";

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  mgr->AddTask(task);
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  
  return task ;
}

