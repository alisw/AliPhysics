AliAnalysisTaskSE * AddTaskPHOSQA()
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C"); 

  AliAnalysisTaskCaloCellsQA *taskPHOSCellQA = AddTaskCaloCellsQA(5, 1, "CaloCellsQA.root","PHOSCellsQA");
  taskPHOSCellQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0); 

  return taskPHOSCellQA;
}
