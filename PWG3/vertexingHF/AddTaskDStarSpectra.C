//DEFINITION OF A FEW CONSTANTS
const Int_t    minITSClusters = 5;
const Int_t    minITSClustersSoft = 4;
const Int_t    numberOfSigmasPID = 3;
// ANALYSIS TYPE DATA/MC
const Bool_t usePIDforKaons = kFALSE;
//----------------------------------------------------

AliAnalysisTaskSEDStarSpectra *AddTaskDStarSpectra(Bool_t theMCon=kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDStarJets", "No analysis manager to connect to.");
    return NULL;
  }  
  
  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliAnalysisTaskSEDStarSpectra *task = new AliAnalysisTaskSEDStarSpectra("AliAnalysisTaskSEDStarSpectra");
  task->SetMinITSClusters(minITSClusters);
  task->SetMinITSClustersSoft(minITSClustersSoft);
  task->SetMinITSClustersSoft(numberOfSigmasPID);
  task->SetMC(theMCon);
  task->SetMC(usePIDforKaons);

  // Create and connect containers for input/output
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DStarSpectra";

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

