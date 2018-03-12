AliAnalysisTaskSE * AddTaskPhysPHOSQA()
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsPhysQA.C"); 

  AliAnalysisTaskCaloCellsPhysQA *taskPHOSCellQA = AddTaskCaloCellsPhysQA(5, 1, "CaloCellsQA.root","PHOSCellsQA");
  taskPHOSCellQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0); 

  AliCaloCellsPhysQA * qa = (AliCaloCellsPhysQA *) taskPHOSCellQA->GetCaloCellsQA();
  qa->SetPhysicsClusterCut(0.3, 3, 12.5e-9); 

  return taskPHOSCellQA;
}