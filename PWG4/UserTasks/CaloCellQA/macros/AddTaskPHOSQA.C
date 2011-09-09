AliAnalysisTaskSE * AddTaskPHOSQA()
{
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/UserTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C"); 

  AliAnalysisTaskCaloCellsQA *taskPHOSCellQA = AddTaskCaloCellsQA(4, 1, "CaloCellsQA.root","PHOSCellsQA");
  taskPHOSCellQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0); 

  return taskPHOSCellQA;
}
