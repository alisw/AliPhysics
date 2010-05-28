//DEFINITION OF A FEW CONSTANTS
const Int_t    minITSClusters     = 5;
const Int_t    minITSClustersSoft = 4;
const Int_t    numberOfSigmasPID  = 3;
// ANALYSIS TYPE 
const Bool_t anaType   = 0;//0 HD; 1 UU;
const Bool_t usePID = kFALSE;
const Bool_t theMCon=kFALSE;
//----------------------------------------------------

AliAnalysisTaskSEDStarSpectra *AddTaskDStarSpectra()
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDStarSpectra", "No analysis manager to connect to.");
    return NULL;
  }  
  
  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliAnalysisTaskSEDStarSpectra *task = new AliAnalysisTaskSEDStarSpectra("AliAnalysisTaskSEDStarSpectra");
  task->SetAnalysisType(anaType);
  task->SetMinITSClusters(minITSClusters);
  task->SetMinITSClustersSoft(minITSClustersSoft);
  task->SetNSigmasPID(numberOfSigmasPID);
  task->SetMC(theMCon);
  task->SetPID(usePID);
  task->SetDebugLevel(0);

  mgr->AddTask(task);

  // Create and connect containers for input/output
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if (anaType == 0) outputfile += ":PWG3_D2H_DStarSpectraHD";
  if (anaType == 1) outputfile += ":PWG3_D2H_DStarSpectraUU";

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputDStar1 = mgr->CreateContainer("DStarSpectrum",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  AliAnalysisDataContainer *coutputDStar2 = mgr->CreateContainer("DStarAll",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  AliAnalysisDataContainer *coutputDStar3 = mgr->CreateContainer("DStarPID3",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  AliAnalysisDataContainer *coutputDStar4 = mgr->CreateContainer("DStarPID2",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  AliAnalysisDataContainer *coutputDStar5 = mgr->CreateContainer("DStarPID1",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutputDStar1);
  mgr->ConnectOutput(task,3,coutputDStar2);
  mgr->ConnectOutput(task,4,coutputDStar3);
  mgr->ConnectOutput(task,5,coutputDStar4);
  mgr->ConnectOutput(task,6,coutputDStar5);

  
  return task ;
}

