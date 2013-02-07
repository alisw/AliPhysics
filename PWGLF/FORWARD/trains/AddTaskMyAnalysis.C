/**
 * @file   AddTaskMyAnalysis.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Nov 23 02:11:39 2012
 * 
 * @brief  Example script to add a task 
 */
/** 
 * @defgroup pwglf_forward_trains_example TrainSetup Examples 
 * 
 * @ingroup pwglf_forward_trains 
 */
/** 
 * Create and add an analysis task to the train 
 * 
 * @ingroup pwglf_forward_trains_examples
 *
 * @return Created analysis task
 */
AliAnalysisTask* AddTaskMyAnalysis()
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  MyAnalysis* task = new MyAnalysis("test");
  mgr->AddTask(task);
    
  AliAnalysisDataContainer* sums = 
    mgr->CreateContainer("Sums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* results = // Needed for output from Terminate
    mgr->CreateContainer("Results", TList::Class(), 
			 AliAnalysisManager::kParamContainer, // Important!
			 AliAnalysisManager::GetCommonFileName());
  
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, results);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  return task;
}
/*
 * EOF
 */
